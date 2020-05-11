package pipeline

import (
	"fmt"
	"io"
	"os"
	"sync"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/seqio"
	"github.com/will-rowe/groot/src/version"
)

// theBoss is used to orchestrate the minions
type theBoss struct {
	info                *Info                    // the runtime info for the pipeline
	graphMinionRegister []*graphMinion           // used to keep a record of the graph minions
	refSAMheaders       map[int][]*sam.Reference // map of SAM headers for each reference sequence, indexed by path ID
	reads               chan *seqio.FASTQread    // the boss uses this channel to receive data from the main sketching pipeline
	alignments          chan *sam.Record         // used to receive alignments from the graph minions
	bamwriter           *bam.Writer              // destination for the BAM output
	receivedReadCount   int                      // the number of reads the boss is sent during it's lifetime
	mappedCount         int                      // the total number of reads that were successful mapped to at least one graph
	multimappedCount    int                      // the total number of reads that had mappings to multiple graphs
	alignmentCount      int                      // the total number of alignment segments reported post hierarchical alignment of mapped reads
	sync.Mutex                                   // allows sketching minions to update the Boss's count
}

// newBoss will initialise and return theBoss
func newBoss(runtimeInfo *Info, inputChan chan *seqio.FASTQread) *theBoss {
	return &theBoss{
		info:              runtimeInfo,
		reads:             inputChan,
		alignments:        make(chan *sam.Record, BUFFERSIZE),
		receivedReadCount: 0,
		mappedCount:       0,
		multimappedCount:  0,
		alignmentCount:    0,
	}
}

// setupBAM will set up the BAM STDOUT for reporting exact graph alignments
func (theBoss *theBoss) setupBAM() error {

	// get the SAM headers
	samHeaders, err := theBoss.info.Store.GetSAMrefs()
	if err != nil {
		return err
	}
	theBoss.refSAMheaders = samHeaders

	// get program info for SAM header (unique ID, name, command, previous program ID, version)
	programInfo := sam.NewProgram("1", "groot", "groot align", "", version.GetVersion())

	// get some readgroup information TODO: set this properly
	rg, err := sam.NewReadGroup("readsID", "", "", "", "groot align", "illumina", "", "sampleID", "", "", time.Now(), 1000)
	if err != nil {
		return err
	}

	// get all the reference sequences ready for the SAM file
	references := []*sam.Reference{}
	for _, refs := range samHeaders {
		references = append(references, refs...)
	}

	// create the header
	header, err := sam.NewHeader(nil, references)
	if err != nil {
		return err
	}
	header.Version = "1.5"

	// add the program info
	if err := header.AddProgram(programInfo); err != nil {
		return err
	}

	// add the readgroup info
	if err := header.AddReadGroup(rg); err != nil {
		return err
	}

	// use a BAM file or STDOUT (TODO: not exposed to CLI yet)
	var fh io.Writer
	if theBoss.info.Sketch.BAMout != "" {
		var err error
		fh, err = os.Create(theBoss.info.Sketch.BAMout)
		if err != nil {
			return (fmt.Errorf("could not open file for BAM writing: %v", err))
		}
	} else {
		fh = os.Stdout
	}

	// create the bam writer and write the header
	bw, err := bam.NewWriter(fh, header, 0)
	if err != nil {
		return err
	}
	theBoss.bamwriter = bw
	return nil
}

// mapReads is a method to start off the minions to map and align reads, the minions to augment graphs, and collate the alignments
func (theBoss *theBoss) mapReads() error {
	theBoss.alignments = make(chan *sam.Record, BUFFERSIZE)

	// set up the BAM if exact alignment is requested
	if !theBoss.info.Sketch.NoExactAlign {
		if err := theBoss.setupBAM(); err != nil {
			return err
		}
	}

	// setup the waitgroups for the sketching and graphing minions
	var wg1 sync.WaitGroup
	var wg2 sync.WaitGroup

	// launch the graph minions (one minion per graph in the index)
	theBoss.graphMinionRegister = make([]*graphMinion, len(theBoss.info.Store))
	for _, graph := range theBoss.info.Store {

		// create, start and register the graph minion
		minion := newGraphMinion(theBoss, graph)
		wg2.Add(1)
		minion.start(&wg2)
		theBoss.graphMinionRegister[graph.GraphID] = minion
	}

	// launch the sketching minions (one per CPU)
	for i := 0; i < theBoss.info.NumProc; i++ {
		wg1.Add(1)
		go func(workerNum int) {
			defer wg1.Done()

			// keep a track of what this minion does
			receivedReads := 0
			mappedCount := 0
			multimappedCount := 0

			// start the main processing loop
			for {

				// pull reads from queue until done
				read, ok := <-theBoss.reads
				if !ok {

					// update the counts
					theBoss.Lock()
					theBoss.receivedReadCount += receivedReads
					theBoss.mappedCount += mappedCount
					theBoss.multimappedCount += multimappedCount
					theBoss.Unlock()

					// end the sketching minion
					return
				}

				// get sketch for read
				readSketch, err := read.RunMinHash(theBoss.info.KmerSize, theBoss.info.SketchSize, false, nil)
				if err != nil {
					panic(err)
				}

				// get the number of k-mers in the sequence
				kmerCount := (len(read.Seq) - theBoss.info.KmerSize) + 1

				// query the LSH ensemble
				results, err := theBoss.info.db.Query(readSketch, kmerCount, theBoss.info.ContainmentThreshold)
				if err != nil {
					panic(err)
				}

				// if multiple graphs are returned, we need to deep copy the read
				deepCopy := false
				if len(results) > 1 {
					deepCopy = true
				}

				// augment graphs and optionally perform exact alignment
				for graphID, hits := range results {
					if deepCopy {
						readCopy := *read.DeepCopy()
						theBoss.graphMinionRegister[graphID].inputChannel <- &graphMinionPair{hits, readCopy}
					} else {
						theBoss.graphMinionRegister[graphID].inputChannel <- &graphMinionPair{hits, *read}
					}
				}

				// update counts
				receivedReads++
				if len(results) > 0 {
					mappedCount++
				}
				if len(results) > 1 {
					multimappedCount++
				}
			}
		}(i)
	}

	// control the channels
	go func() {

		// wait for the sketching minions to finish
		wg1.Wait()

		// shut down the graph minions input channels
		for _, minion := range theBoss.graphMinionRegister {
			close(minion.inputChannel)
		}

		// wait for the graph minions to finish
		wg2.Wait()

		// end the alignment writer
		close(theBoss.alignments)

	}()

	// collect the alignments and write them
	for record := range theBoss.alignments {
		// check the record is valid
		//if sam.IsValidRecord(record) == false {
		//	os.Exit(1)
		//}
		theBoss.alignmentCount++
		if err := theBoss.bamwriter.Write(record); err != nil {
			return err
		}
	}

	// close the bam writer and return to the completed boss to the pipeline
	var err error
	if !theBoss.info.Sketch.NoExactAlign {
		err = theBoss.bamwriter.Close()
	}
	return err
}
