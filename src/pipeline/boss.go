package pipeline

import (
	"log"
	"os"
	"sync"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/lshforest"
	"github.com/will-rowe/groot/src/minhash"
	"github.com/will-rowe/groot/src/misc"
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
	receivedReadCount   int                      // the number of reads the boss is sent during it's lifetime
	mappedCount         int                      // the total number of reads that were successful mapped to at least one graph
	multimappedCount    int                      // the total number of reads that had multiple mappings
	alignmentCount      int                      // the total number of alignment segments reported post hierarchical alignment of mapped reads
	sync.Mutex                                   // allows sketching minions to update the Boss's count
}

// newBoss will initialise and return theBoss
func newBoss(runtimeInfo *Info, inputChan chan *seqio.FASTQread) *theBoss {
	return &theBoss{
		info:              runtimeInfo,
		reads:             inputChan,
		receivedReadCount: 0,
		mappedCount:       0,
		multimappedCount:  0,
		alignmentCount:    0,
	}
}

// mapReads is a method to start off the minions to map and align reads, the minions to augment graphs, and collate the alignments
func (theBoss *theBoss) mapReads() error {

	// get the SAM headers
	samHeaders, err := theBoss.info.Store.GetSAMrefs()
	if err != nil {
		return err
	}
	theBoss.refSAMheaders = samHeaders
	theBoss.alignments = make(chan *sam.Record, BUFFERSIZE)

	// get program info for SAM header (unique ID, name, command, previous program ID, version)
	programInfo := sam.NewProgram("1", "groot", "groot align", "", version.VERSION)

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

	// create the bam writer and write the header
	bw, err := bam.NewWriter(os.Stdout, header, 0)
	if err != nil {
		return err
	}

	// setup the waitgroups for the sketching and graphing minions
	var wg1 sync.WaitGroup
	var wg2 sync.WaitGroup

	// launch the graph minions (one minion per graph in the index)
	theBoss.graphMinionRegister = make([]*graphMinion, len(theBoss.info.Store))
	for graphID, graph := range theBoss.info.Store {

		// create, start and register the graph minion
		minion := newGraphMinion(graphID, graph, theBoss.alignments, theBoss.refSAMheaders[int(graphID)])
		wg2.Add(1)
		minion.start(&wg2)
		theBoss.graphMinionRegister[graphID] = minion
	}
	log.Printf("\tgraph workers: %d", len(theBoss.graphMinionRegister))

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
				readSketch, err := minhash.GetReadSketch(read.Seq, uint(theBoss.info.KmerSize), uint(theBoss.info.SketchSize), false)
				misc.ErrorCheck(err)

				// get the number of k-mers in the sequence
				kmerCount := (len(read.Seq) - theBoss.info.KmerSize) + 1

				// query the LSH ensemble
				hits, err := theBoss.info.db.Query(readSketch, kmerCount, theBoss.info.ContainmentThreshold)
				if err != nil {
					panic(err)
				}
				for _, hit := range hits {

					// make a copy of this graphWindow
					graphWindow := &lshforest.Key{
						GraphID:        hit.GraphID,
						Node:           hit.Node,
						OffSet:         hit.OffSet,
						ContainedNodes: hit.ContainedNodes, // don't need to deep copy this as we don't edit it
						Freq:           float64(kmerCount), // add the k-mer count of the read in this window
					}

					// send the window to the correct go routine for read alignment and graph augmentation
					theBoss.graphMinionRegister[hit.GraphID].inputChannel <- graphMinionPair{graphWindow, read}
				}

				// update counts
				receivedReads++
				if len(hits) > 0 {
					mappedCount++
				}
				if len(hits) > 1 {
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
		if err := bw.Write(record); err != nil {
			return err
		}
	}

	// close the bam writer and return to the completed boss to the pipeline
	return bw.Close()
}
