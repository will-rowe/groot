package pipeline

/*
 this part of the pipeline will process reads, sketch them, map them and then project them onto variation graphs
*/

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"log"
	"os"
	"strings"
	"sync"

	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
)

// DataStreamer is a pipeline process that streams data from STDIN/file
type DataStreamer struct {
	info   *Info
	input  []string
	output chan []byte
}

// NewDataStreamer is the constructor
func NewDataStreamer(info *Info) *DataStreamer {
	return &DataStreamer{info: info, output: make(chan []byte, BUFFERSIZE)}
}

// Connect is the method to connect the DataStreamer to some data source
func (proc *DataStreamer) Connect(input []string) {
	proc.input = input
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *DataStreamer) Run() {
	defer close(proc.output)
	var scanner *bufio.Scanner
	// if an input file path has not been provided, scan the contents of STDIN
	if len(proc.input) == 0 {
		scanner = bufio.NewScanner(os.Stdin)
		for scanner.Scan() {
			// important: copy content of scan to a new slice before sending, this avoids race conditions (as we are using multiple go routines) from concurrent slice access
			proc.output <- append([]byte(nil), scanner.Bytes()...)
		}
		if scanner.Err() != nil {
			log.Fatal(scanner.Err())
		}
	} else {
		for i := 0; i < len(proc.input); i++ {
			fh, err := os.Open(proc.input[i])
			misc.ErrorCheck(err)
			defer fh.Close()
			// handle gzipped input
			splitFilename := strings.Split(proc.input[i], ".")
			if splitFilename[len(splitFilename)-1] == "gz" {
				gz, err := gzip.NewReader(fh)
				misc.ErrorCheck(err)
				defer gz.Close()
				scanner = bufio.NewScanner(gz)
			} else {
				scanner = bufio.NewScanner(fh)
			}
			for scanner.Scan() {
				proc.output <- append([]byte(nil), scanner.Bytes()...)
			}
			if scanner.Err() != nil {
				log.Fatal(scanner.Err())
			}
		}
	}
}

// WASMstreamer is a pipeline process that streams data from the WASM JS function
type WASMstreamer struct {
	input  chan []byte
	output chan []byte
}

// NewWASMstreamer is the constructor
func NewWASMstreamer() *WASMstreamer {
	return &WASMstreamer{output: make(chan []byte, BUFFERSIZE)}
}

// ConnectChan is a to connect the pipeline to the WASM JS function
func (proc *WASMstreamer) ConnectChan(inputChan chan []byte) {
	proc.input = inputChan
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *WASMstreamer) Run() {
	defer close(proc.output)
	leftOvers := []byte{}

	// collect a chunk of fastq from the WASM JS function
	for chunk := range proc.input {

		// ignore any empty receives
		if len(chunk) == 0 {
			continue
		}

		// add any leftovers from the last chunk to the start of this chunk
		if len(leftOvers) != 0 {
			chunk = append(leftOvers, chunk...)
			leftOvers = nil
		}

		// check the current chunk to see if it ends with a truncated line
		// if it does, find the final line break and create subchunk
		i := len(chunk) - 1
		if chunk[i] != 0x0A {
			for i > 0 {
				i--
				if chunk[i] == 0x0A {
					break
				}
			}
			leftOvers = chunk[i:]
			if i == 0 {
				continue
			}
			chunk = chunk[0:i]
		}

		// convert bytes to reader
		chunkBuffer := bytes.NewReader(chunk)

		// read the chunk line by line
		scanner := bufio.NewScanner(chunkBuffer)
		for scanner.Scan() {

			// ignore empty lines and send
			line := bytes.TrimSpace(scanner.Bytes())

			if len(line) > 0 {
				proc.output <- append([]byte(nil), line...)
			}
		}
		if scanner.Err() != nil {
			log.Fatal(scanner.Err())
		}

	}
}

// FastqHandler is a pipeline process to convert a pipeline to the FASTQ type
type FastqHandler struct {
	info   *Info
	input  chan []byte
	output chan *seqio.FASTQread
}

// NewFastqHandler is the constructor
func NewFastqHandler(info *Info) *FastqHandler {
	return &FastqHandler{info: info, output: make(chan *seqio.FASTQread, BUFFERSIZE)}
}

// ConnectWASM  is a tmp solution for WASM
func (proc *FastqHandler) ConnectWASM(previous *WASMstreamer) {
	proc.input = previous.output
}

// Connect is the method to join the input of this process with the output of a DataStreamer
func (proc *FastqHandler) Connect(previous *DataStreamer) {
	proc.input = previous.output
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *FastqHandler) Run() {
	defer close(proc.output)
	var l1, l2, l3, l4 []byte
	if proc.info.Sketch.Fasta {
		for line := range proc.input {
			if len(line) == 0 {
				break
			}

			// check for chevron
			if line[0] == 62 {
				if l1 != nil {

					// store current fasta entry (as FASTQ read)
					l1[0] = 64
					newRead, err := seqio.NewFASTQread(l1, l2, nil, nil)
					if err != nil {
						log.Fatal(err)
					}

					// send on the new read and reset the line stores
					proc.output <- newRead
				}
				l1, l2 = line, nil
			} else {
				l2 = append(l2, line...)
			}
		}

		// flush final fasta
		l1[0] = 64
		newRead, err := seqio.NewFASTQread(l1, l2, nil, nil)
		if err != nil {
			log.Fatal(err)
		}

		// send on the new read and reset the line stores
		proc.output <- newRead
	} else {

		// grab four lines and create a new FASTQread struct from them - perform some format checks and trim low quality bases
		for line := range proc.input {
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
				proc.output <- newRead
				l1, l2, l3, l4 = nil, nil, nil, nil
			}
		}
	}
}

// FastqChecker is a process to quality check FASTQ reads and send the sequence on for mapping
type FastqChecker struct {
	info   *Info
	input  chan *seqio.FASTQread
	output chan *seqio.FASTQread
}

// NewFastqChecker is the constructor
func NewFastqChecker(info *Info) *FastqChecker {
	return &FastqChecker{info: info, output: make(chan *seqio.FASTQread, BUFFERSIZE)}
}

// Connect is the method to join the input of this process with the output of FastqHandler
func (proc *FastqChecker) Connect(previous *FastqHandler) {
	proc.input = previous.output
}

// Run is the method to run this process, which satisfies the pipeline interface
// TODO: I've removed the QC bits for now
func (proc *FastqChecker) Run() {
	log.Printf("now streaming reads...")

	// count the number of reads and their lengths as we go
	rawCount, lengthTotal := 0, 0
	for read := range proc.input {
		rawCount++

		// tally the length so we can report the mean
		lengthTotal += len(read.Seq)

		// send the read onwards for mapping
		proc.output <- read
	}

	// check we have received reads & print stats
	if rawCount == 0 {
		misc.ErrorCheck(errors.New("no fastq reads received"))
	}
	log.Printf("\tnumber of reads received from input: %d\n", rawCount)
	meanRL := float64(lengthTotal) / float64(rawCount)
	log.Printf("\tmean read length: %.0f\n", meanRL)
	close(proc.output)
}

// ReadMapper is a pipeline process to query the LSH database, map reads and project alignments onto graphs
type ReadMapper struct {
	info      *Info
	input     chan *seqio.FASTQread
	output    chan *graph.GrootGraph
	readStats [4]int // corresponds to num. reads, total num. mapped, num. multimapped, total k-mers
}

// NewReadMapper is the constructor
func NewReadMapper(info *Info) *ReadMapper {
	return &ReadMapper{info: info, output: make(chan *graph.GrootGraph), readStats: [4]int{0, 0, 0, 0}}
}

// Connect is the method to join the input of this process with the output of FastqChecker
func (proc *ReadMapper) Connect(previous *FastqChecker) {
	proc.input = previous.output
}

// CollectReadStats is a method to return the number of reads processed, how many mapped and the number of multimaps
func (proc *ReadMapper) CollectReadStats() [4]int {
	return proc.readStats
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *ReadMapper) Run() {
	defer close(proc.output)

	// set up the boss/minion pool and run the mapping
	theBoss := newBoss(proc.info, proc.input)
	misc.ErrorCheck(theBoss.mapReads())

	// log some stuff
	if theBoss.receivedReadCount == 0 {
		misc.ErrorCheck(fmt.Errorf("no reads passed quality-based trimming"))
	} else {
		log.Printf("\tnumber of reads sketched: %d\n", theBoss.receivedReadCount)
	}
	proc.readStats[0] = theBoss.receivedReadCount
	proc.readStats[1] = theBoss.mappedCount
	proc.readStats[2] = theBoss.multimappedCount

	// nothing may have mapped, which isn't an error - so make GROOT exit gracefully
	if proc.readStats[1] == 0 {
		log.Println("no reads could be mapped to the reference graphs")

		// set the graph store to nil so we don't end up writing anything
		proc.info.Store = make(graph.Store)
		return
	}
	log.Printf("\ttotal number of mapped reads: %d\n", theBoss.mappedCount)
	log.Printf("\t\tuniquely mapped: %d\n", (theBoss.mappedCount - theBoss.multimappedCount))
	log.Printf("\t\tmultimapped: %d\n", theBoss.multimappedCount)
	log.Printf("\ttotal number of alignments: %d\n", theBoss.alignmentCount)

	// send on the graphs for pruning now that the mapping is done
	for _, g := range proc.info.Store {
		proc.readStats[3] += int(g.KmerTotal)
		proc.output <- g
	}
	log.Print("processing graphs...")
	log.Printf("\ttotal number of k-mers projected onto graphs: %d\n", proc.readStats[3])

	// fix for Wasm: needs total number of processed kmers recording as we don't recalc this later
	proc.info.Haplotype.TotalKmers = proc.readStats[3]
}

// GraphPruner is a pipeline process to prune the graphs post mapping
type GraphPruner struct {
	info             *Info
	input            chan *graph.GrootGraph
	output           chan *graph.GrootGraph
	foundPaths       []string
	connectHaplotype bool // activate this flag if this pipeline process is connected to the EMpathFinder
}

// NewGraphPruner is the constructor
func NewGraphPruner(info *Info, conH bool) *GraphPruner {
	return &GraphPruner{info: info, output: make(chan *graph.GrootGraph, BUFFERSIZE), connectHaplotype: conH}
}

// Connect is the method to join the input of this process with the output of ReadMapper
func (proc *GraphPruner) Connect(previous *ReadMapper) {
	proc.input = previous.output
}

// CollectOutput is a method to return what paths are left post-pruning
func (proc *GraphPruner) CollectOutput() []string {
	return proc.foundPaths
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *GraphPruner) Run() {
	defer close(proc.output)
	graphChan := make(chan *graph.GrootGraph)
	var wg sync.WaitGroup
	counter := 0
	for g := range proc.input {
		wg.Add(1)
		counter++

		go func(graph *graph.GrootGraph) {
			defer wg.Done()
			// check for alignments and prune the graph
			keepGraph := graph.Prune(proc.info.Sketch.MinKmerCoverage)

			// check we have some graph
			if keepGraph != false {
				graphChan <- graph
			}
		}(g)
	}
	go func() {
		wg.Wait()
		close(graphChan)
	}()

	// count and print some stuff
	keptPaths := []string{}
	keptGraphs := make(graph.Store)
	for g := range graphChan {
		g.GrootVersion = proc.info.Version
		keptGraphs[g.GraphID] = g
		log.Printf("\tgraph %d has %d remaining paths after weighting and pruning", g.GraphID, len(g.Paths))
		for _, path := range g.Paths {
			log.Printf("\t- [%v]", string(path))
			keptPaths = append(keptPaths, string(path))
		}
		if proc.connectHaplotype {
			proc.output <- g
		}
	}
	if counter == 0 {
		return
	}
	log.Printf("\ttotal number of graphs pruned: %d\n", counter)
	if len(keptGraphs) == 0 {
		log.Print("\tno graphs remaining after pruning")
		return
	}
	log.Printf("\ttotal number of graphs remaining: %d\n", len(keptGraphs))
	log.Printf("\ttotal number of possible haplotypes found: %d\n", len(keptPaths))
	proc.info.Store = keptGraphs
	proc.foundPaths = keptPaths
}
