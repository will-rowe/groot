package pipeline

import (
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/lshforest"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
)

type graphMinionPair struct {
	mapping *lshforest.Key
	read    *seqio.FASTQread
}

// graphMinion holds a graph and is responsible for augmenting the paths when new mapping data arrives
type graphMinion struct {
	id            uint32
	graph         *graph.GrootGraph
	inputChannel  chan graphMinionPair
	outputChannel chan *sam.Record
	runAlignment  bool
	references    []*sam.Reference // the SAM references for each path in this graph
}

// newGraphMinion is the constructor function
func newGraphMinion(id uint32, graph *graph.GrootGraph, alignmentChan chan *sam.Record, refs []*sam.Reference) *graphMinion {
	return &graphMinion{
		id:            id,
		graph:         graph,
		inputChannel:  make(chan graphMinionPair, BUFFERSIZE),
		outputChannel: alignmentChan,
		runAlignment:  true, // setting this to True for now to replicate groot behaviour - can be used later to run groot haplotype workflow etc.
		references:    refs,
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

			// increment the graph node weightings for nodes contained in the mapping window
			misc.ErrorCheck(graphMinion.graph.IncrementSubPath(mappingData.mapping.ContainedNodes, mappingData.mapping.Freq))

			// perform the alignment if requested
			if !graphMinion.runAlignment {
				continue
			}

			// perform graph alignment on forward and reverse complement
			// TODO: as we used canonical k-mer hashing, not sure which orientation the read seeded in, think about this some more
			for i := 0; i < 2; i++ {

				// run the alignment
				alignments, err := graphMinion.graph.AlignRead(mappingData.read, mappingData.mapping, graphMinion.references)
				if err != nil {
					panic(err)
				}
				for _, alignment := range alignments {
					graphMinion.outputChannel <- alignment
				}

				// reverse complement read and run again
				mappingData.read.RevComplement()
				if i == 1 {
					break
				}
			}
		}
	}()
}
