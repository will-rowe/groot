package pipeline

import (
	"sort"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/lshe"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
)

type graphMinionPair struct {
	mappings lshe.Keys
	read     seqio.FASTQread
}

// graphMinion holds a graph and is responsible for augmenting the paths when new mapping data arrives
type graphMinion struct {
	id            uint32
	graph         *graph.GrootGraph
	inputChannel  chan *graphMinionPair
	outputChannel chan *sam.Record
	runAlignment  bool
	references    []*sam.Reference // the SAM references for each path in this graph
	boss          *theBoss         // pointer to the boss so the minion can access the runtime info (e.g. k-mer size)
}

// newGraphMinion is the constructor function
func newGraphMinion(id uint32, graph *graph.GrootGraph, alignmentChan chan *sam.Record, refs []*sam.Reference, boss *theBoss) *graphMinion {
	return &graphMinion{
		id:            id,
		graph:         graph,
		inputChannel:  make(chan *graphMinionPair, BUFFERSIZE),
		outputChannel: alignmentChan,
		references:    refs,
		boss:          boss,
	}
}

// start is a method to start the graphMinion running
func (graphMinion *graphMinion) start(wg *sync.WaitGroup) {
	go func() {
		for {

			// pull read mappings from queue until done
			mappingData, ok := <-graphMinion.inputChannel
			if !ok {
				wg.Done()
				return
			}

			// sort the mappings for this read
			sort.Sort(mappingData.mappings)

			// claculate the number of k-mers for the read
			kmerCount := float64(len(mappingData.read.Seq)-graphMinion.boss.info.KmerSize) + 1

			// process each mapping until an exact alignment found
			alignmentFound := false
			for _, mapping := range mappingData.mappings {

				// increment the graph node weightings for nodes contained in the mapping window
				misc.ErrorCheck(graphMinion.graph.IncrementSubPath(mapping.ContainedNodes, kmerCount))

				// perform the alignment if requested
				if graphMinion.boss.info.Sketch.NoExactAlign {
					continue
				}

				// perform graph alignment on forward and reverse complement
				// TODO: as we used canonical k-mer hashing, not sure which orientation the read seeded in, think about this some more
				for i := 0; i < 2; i++ {

					// run the alignment
					alignments, err := graphMinion.graph.AlignRead(&mappingData.read, &mapping, graphMinion.references)
					if err != nil {
						panic(err)
					}

					// if an alignment was found, send them and call it a day
					if len(alignments) != 0 {
						for _, alignment := range alignments {
							graphMinion.outputChannel <- alignment
						}
						alignmentFound = true
						break
					}

					// reverse complement read and run again if no alignment found
					mappingData.read.RevComplement()
				}
				if alignmentFound {
					break
				}
			}
			//log.Printf("graph %d could not find alignment for %v after trying %d mapping locations", graphMinion.id, string(mappingData.read.ID), len(mappingData.mappings))
		}
	}()
}
