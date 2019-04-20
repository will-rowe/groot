/*
	the stream package contains a streaming implementation based on the Gopher Academy article by S. Lampa - Patterns for composable concurrent pipelines in Go (https://blog.gopheracademy.com/advent-2015/composable-pipelines-improvements/)
*/
package stream

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"log"
	"os"
	"strings"
	"sync"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/alignment"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/lshForest"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
	"github.com/will-rowe/groot/src/version"
)

const (
	BUFFERSIZE = 128 // buffer size to use for channels
)

/*
  The process interface
*/
type process interface {
	Run()
}

/*
  The basic pipeline - takes a list of Processes and runs them in Go routines, the last process is ran in the fg
*/
type Pipeline struct {
	Processes []process
}

func NewPipeline() *Pipeline {
	return &Pipeline{}
}

func (pl *Pipeline) AddProcess(proc process) {
	pl.Processes = append(pl.Processes, proc)
}

func (pl *Pipeline) AddProcesses(procs ...process) {
	for _, proc := range procs {
		pl.AddProcess(proc)
	}
}

func (pl *Pipeline) Run() {
	for i, proc := range pl.Processes {
		if i < len(pl.Processes)-1 {
			go proc.Run()
		} else {
			proc.Run()
		}
	}
}

/*
  A process to stream data from STDIN/file
*/
type DataStreamer struct {
	process
	Output    chan []byte
	InputFile []string
}

func NewDataStreamer() *DataStreamer {
	return &DataStreamer{Output: make(chan []byte, BUFFERSIZE)}
}

func (proc *DataStreamer) Run() {
	var scanner *bufio.Scanner
	// if an input file path has not been provided, scan the contents of STDIN
	if len(proc.InputFile) == 0 {
		scanner = bufio.NewScanner(os.Stdin)
		for scanner.Scan() {
			// important: copy content of scan to a new slice before sending, this avoids race conditions (as we are using multiple go routines) from concurrent slice access
			proc.Output <- append([]byte(nil), scanner.Bytes()...)
		}
		if scanner.Err() != nil {
			log.Fatal(scanner.Err())
		}
	} else {
		for i := 0; i < len(proc.InputFile); i++ {
			fh, err := os.Open(proc.InputFile[i])
			misc.ErrorCheck(err)
			defer fh.Close()
			// handle gzipped input
			splitFilename := strings.Split(proc.InputFile[i], ".")
			if splitFilename[len(splitFilename)-1] == "gz" {
				gz, err := gzip.NewReader(fh)
				misc.ErrorCheck(err)
				defer gz.Close()
				scanner = bufio.NewScanner(gz)
			} else {
				scanner = bufio.NewScanner(fh)
			}
			for scanner.Scan() {
				proc.Output <- append([]byte(nil), scanner.Bytes()...)
			}
			if scanner.Err() != nil {
				log.Fatal(scanner.Err())
			}
		}
	}
	close(proc.Output)
}

/*
  A process to generate a FASTQ read from a stream of bytes
*/
type FastqHandler struct {
	process
	Input  chan []byte
	Output chan seqio.FASTQread
}

func NewFastqHandler() *FastqHandler {
	return &FastqHandler{Output: make(chan seqio.FASTQread, BUFFERSIZE)}
}

func (proc *FastqHandler) Run() {
	defer close(proc.Output)
	var l1, l2, l3, l4 []byte
	// grab four lines and create a new FASTQread struct from them - perform some format checks and trim low quality bases
	for line := range proc.Input {
		if l1 == nil {
			l1 = line
		} else if l2 == nil {
			l2 = line
		} else if l3 == nil {
			l3 = line
		} else if l4 == nil {
			l4 = line
			// create fastq read
			newRead, err := seqio.NewFASTQread(l1, l2, l3, l4)
			if err != nil {
				log.Fatal(err)
			}
			// send on the new read and reset the line stores
			proc.Output <- newRead
			l1, l2, l3, l4 = nil, nil, nil, nil
		}
	}
}

/*
  A process to quality check FASTQ reads - trimming and discarding them according to user supplied cut offs
*/
type FastqChecker struct {
	process
	Input         chan seqio.FASTQread
	Output        chan seqio.FASTQread
	WindowSize    int
	MinReadLength int
	MinQual       int
}

func NewFastqChecker() *FastqChecker {
	return &FastqChecker{Output: make(chan seqio.FASTQread, BUFFERSIZE)}
}

func (proc *FastqChecker) Run() {
	defer close(proc.Output)
	log.Printf("now streaming reads...")
	var wg sync.WaitGroup
	// count the number of reads and their lengths as we go
	rawCount, lengthTotal := 0, 0
	for read := range proc.Input {
		rawCount++
		//  tally the length so we can report the mean
		lengthTotal += len(read.Seq)
		// quality-based trim the read if requested -- moving this in to the alignment routine for now...
		/*if proc.MinQual != 0 {
			wg.Add(1)
			go func(read seqio.FASTQread) {
				defer wg.Done()
				read.QualTrim(proc.MinQual)
				// only send read on if it is greater than the length cutoff
				if len(read.Seq) > proc.MinReadLength {
					proc.Output <- read
				}
			}(read)
			// if trimming wasn't requested, only send read on if it is greater than the length cutoff
		} else {
			proc.Output <- read
		} */
		proc.Output <- read
	}
	wg.Wait()
	// check we have received reads & print stats
	if rawCount == 0 {
		misc.ErrorCheck(errors.New("no fastq reads received"))
	}
	log.Printf("\tnumber of reads received from input: %d\n", rawCount)
	meanRL := float64(lengthTotal) / float64(rawCount)
	log.Printf("\tmean read length: %.0f\n", meanRL)
	// check the length is within +/-10 bases of the graph window
	if meanRL < float64(proc.WindowSize-10) || meanRL > float64(proc.WindowSize+10) {
		misc.ErrorCheck(fmt.Errorf("mean read length is outside the graph window size (+/- 10 bases)\n"))
	}
}

/*
  A process to query the LSH database, perform full MinHash comparisons on top hits and returns putative graph mapping locations
*/
type DbQuerier struct {
	process
	Input       chan seqio.FASTQread
	Output      chan seqio.FASTQread
	Db          *lshForest.LSHforest
	CommandInfo *misc.IndexInfo
	GraphStore  graph.GraphStore
}

func NewDbQuerier() *DbQuerier {
	return &DbQuerier{Output: make(chan seqio.FASTQread, BUFFERSIZE)}
}

func (proc *DbQuerier) Run() {
	defer close(proc.Output)
	// record the number of reads processed by the DbQuerier
	readTally, seededTally := 0, 0
	var wg sync.WaitGroup
	seedCollectionChan := make(chan seqio.FASTQread)
	for read := range proc.Input {
		wg.Add(1)
		go func(read seqio.FASTQread) {
			defer wg.Done()
			seeds := []seqio.Key{}
			// query the read twice
			for rcIterator := 0; rcIterator < 2; rcIterator++ {
				// if this is the second loop through the loop, reverse complement the read
				if rcIterator == 1 {
					read.RevComplement()
					read.RC = true // set RC flag so we can tell which orientation the read is in
				}
				// get signature for read
				readSketch, err := read.RunMinHash(proc.CommandInfo.Ksize, proc.CommandInfo.SigSize)
				misc.ErrorCheck(err)
				// query the LSH forest
				for _, result := range proc.Db.Query(readSketch) {
					seed := proc.Db.KeyLookup[result]
					seed.RC = read.RC
					seeds = append(seeds, seed)
				}
			}
			// send any seeds on to the collector
			if len(seeds) > 0 {
				read.Seeds = seeds
				seedCollectionChan <- read
			}
		}(read)
		readTally++
	}
	// close the channel once all the queries are done
	go func() {
		wg.Wait()
		close(seedCollectionChan)
	}()
	// collect the mapping positions for each read and send the seeded read on for alignment
	for seededRead := range seedCollectionChan {
		seededTally++
		proc.Output <- seededRead
	}

	// log some stuff
	if readTally == 0 {
		misc.ErrorCheck(fmt.Errorf("no reads passed quality-based trimming"))
	} else {
		log.Printf("\tnumber of reads received for alignment post QC: %d\n", readTally)
	}
	if seededTally == 0 {
		misc.ErrorCheck(fmt.Errorf("no reads could be seeded against the reference graphs"))
	} else {
		log.Printf("\tnumber of reads seeded: %d\n", seededTally)
	}
}

/*
  A process to align reads once they have been seeded against a graph
*/
type Aligner struct {
	process
	Input         chan seqio.FASTQread
	Output        chan *sam.Record
	GraphStore    graph.GraphStore
	RefMap        map[int][]*sam.Reference
	MinReadLength int
	MinQual       int
	MaxClip       int
}

func NewAligner() *Aligner {
	return &Aligner{Output: make(chan *sam.Record, BUFFERSIZE)}
}

func (proc *Aligner) Run() {
	defer close(proc.Output)
	alignmentCollectionChan := make(chan *sam.Record)
	var wg sync.WaitGroup
	multiSeeds := 0
	// receive seeded reads
	for seededRead := range proc.Input {
		// check the number of seeds found for the read (i.e. find number of multimappers)
		if len(seededRead.Seeds) > 1 {
			multiSeeds++
		}
		// begin the local alignment
		wg.Add(1)
		go func(read seqio.FASTQread) {
			defer wg.Done()
			readAlignmentCount := 0
			// seededGraphs records each graph seeded and if a local alignment is achieved - it prevents duplicate alignments being reported from closely separated seeds
			seededGraphs := make(map[int]int)
			// quality trim the read if requested
			if proc.MinQual != 0 {
				read.QualTrim(proc.MinQual)
				// only send read on if it is greater than the length cutoff
				if len(read.Seq) < proc.MinReadLength {
					return
				}
			}
			// try alignment on all the seeds
			for i := 0; i < len(read.Seeds); i++ {
				// if the seed required RC of read, make sure the read is in the right complement orientation
				if read.RC != read.Seeds[i].RC {
					read.RevComplement()
				}
				// perform hierarchical local alignment of seeded read against the variation graph
				alignments := alignment.Align(read, i, proc.GraphStore[read.Seeds[i].GraphID], proc.RefMap[read.Seeds[i].GraphID], proc.MaxClip)
				// range over the received alignments channel and send the alignments from this read on to the collection channel
				for alignment := range alignments {
					readAlignmentCount++

					if readAlignmentCount > 1 && read.RC == true {
						alignment.Flags = sam.Secondary | sam.Reverse
					} else if read.RC == true {
						alignment.Flags = sam.Reverse
					} else if readAlignmentCount > 1 {
						alignment.Flags = sam.Secondary
					}

					alignmentCollectionChan <- alignment
				}
				if readAlignmentCount != 0 {
					seededGraphs[read.Seeds[i].GraphID] = read.Seeds[i].Node
				}
			}
			// if no alignments, record read as unmapped (allows us to see what seeded reads couldn't be mapped)
			if readAlignmentCount == 0 {
				samReference, _ := sam.NewReference("*", "", "", 0, nil, nil)
				record := &sam.Record{
					Name:  string(read.ID[1:]),
					Flags: sam.Unmapped,
					Ref:   samReference,
					Pos:   -1,
				}
				// Unmapped reads should be stored in the orientation in which they came off the sequencing machine
				if read.RC == true {
					read.RevComplement()
				}
				record.Seq = sam.NewSeq(read.Seq)
				record.Qual = read.Qual
				alignmentCollectionChan <- record
			}
		}(seededRead)
	}
	go func() {
		wg.Wait()
		close(alignmentCollectionChan)
	}()

	// collect the alignments and send to the SAM writer
	for alignment := range alignmentCollectionChan {
		proc.Output <- alignment
	}
	// report some stuff TODO: add more information here...
	log.Printf("\tnumber of reads with multiple seeds: %d\n", multiSeeds)
}

/*
  A process to drive the completion of the pipeline and output the SAM
*/
type SamWriter struct {
	process
	Input  chan *sam.Record
	RefMap map[int][]*sam.Reference
}

func NewSamWriter() *SamWriter {
	return &SamWriter{}
}

func (proc *SamWriter) Run() {
	// get program info for SAM header (unique ID, name, command, previous program ID, version)
	programInfo := sam.NewProgram("1", "groot", "groot align", "", version.VERSION)
	// get some readgroup information TODO: set this properly
	rg, err := sam.NewReadGroup("readsID", "", "", "", "groot align", "illumina", "", "sampleID", "", "", time.Now(), 1000)
	misc.ErrorCheck(err)
	// get all the reference sequences ready for the SAM file
	references := []*sam.Reference{}
	for _, refs := range proc.RefMap {
		references = append(references, refs...)
	}
	// create the header
	header, err := sam.NewHeader(nil, references)
	misc.ErrorCheck(err)
	header.Version = "1.5"
	// add the program info
	misc.ErrorCheck(header.AddProgram(programInfo))
	// add the readgroup info
	misc.ErrorCheck(header.AddReadGroup(rg))
	// create the bam writer and write the header
	bw, err := bam.NewWriter(os.Stdout, header, 0)
	misc.ErrorCheck(err)
	// collect the alignments and write them
	for record := range proc.Input {
		// check the record is valid
		//if sam.IsValidRecord(record) == false {
		//	os.Exit(1)
		//}
		err := bw.Write(record)
		if err != nil {
			panic(err)
		}
	}
	misc.ErrorCheck(bw.Close())
}
