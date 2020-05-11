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
	boss         *theBoss // pointer to the boss so the minion can access the runtime info (e.g. k-mer size) and channels etc
	id           uint32   // id corresponds to the graphID
	graph        *graph.GrootGraph
	inputChannel chan *graphMinionPair
	runAlignment bool
	references   []*sam.Reference // the SAM references for each path in this graph
}

// newGraphMinion is the constructor function
func newGraphMinion(boss *theBoss, graph *graph.GrootGraph) *graphMinion {
	return &graphMinion{
		boss:         boss,
		id:           graph.GraphID,
		graph:        graph,
		inputChannel: make(chan *graphMinionPair, BUFFERSIZE),
	}
}

// start is a method to start the graphMinion running
func (graphMinion *graphMinion) start(wg *sync.WaitGroup) {

	// get the correct SAM formatted refs for this graph
	graphMinion.references = graphMinion.boss.refSAMheaders[int(graphMinion.id)]

	//
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
							graphMinion.boss.alignments <- alignment
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
