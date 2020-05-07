package pipeline

/*
 this part of the pipeline will convert multiple sequence alignments to variation graphs, sketch the traversals and index them
*/

import (
	"fmt"
	"log"
	"sync"

	"github.com/biogo/biogo/seq/multi"
	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/lshe"
	"github.com/will-rowe/groot/src/misc"
)

// MSAconverter is a pipeline process that converts a list of MSAs to GFAs
type MSAconverter struct {
	info   *Info
	input  []string
	output chan *graph.GrootGraph
}

// NewMSAconverter is the constructor
func NewMSAconverter(info *Info) *MSAconverter {
	return &MSAconverter{info: info, output: make(chan *graph.GrootGraph, BUFFERSIZE)}
}

// Connect is the method to connect the MSAconverter to some data source
func (proc *MSAconverter) Connect(input []string) {
	proc.input = input
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *MSAconverter) Run() {
	var wg sync.WaitGroup
	wg.Add(len(proc.input))

	// load each MSA outside of the go-routines to prevent 'too many open files' error on OSX
	for i, msaFile := range proc.input {
		msa, err := gfa.ReadMSA(msaFile)
		misc.ErrorCheck(err)
		go func(msaID int, msa *multi.Multi) {
			defer wg.Done()

			// convert the MSA to a GFA instance
			newGFA, err := gfa.MSA2GFA(msa)
			misc.ErrorCheck(err)

			// create a GrootGraph
			grootGraph, err := graph.CreateGrootGraph(newGFA, msaID)
			if err != nil {
				misc.ErrorCheck(err)
			}

			// mark the graph has masked if the requested window size is larger than the smallest seq in the graph
			for i, seqLen := range grootGraph.Lengths {
				if seqLen < proc.info.WindowSize {
					log.Printf("\tsequence for %v is shorter than window size (%d vs. %d), skipping graph", string(grootGraph.Paths[i]), seqLen, proc.info.WindowSize)
					grootGraph.Masked = true
					break
				}
			}
			proc.output <- grootGraph
		}(i, msa)
	}
	wg.Wait()
	close(proc.output)
}

// GraphSketcher is a pipeline process that windows graph traversals and sketches them
type GraphSketcher struct {
	info   *Info
	input  chan *graph.GrootGraph
	output chan map[string][]lshe.Key
}

// NewGraphSketcher is the constructor
func NewGraphSketcher(info *Info) *GraphSketcher {
	return &GraphSketcher{info: info, output: make(chan map[string][]lshe.Key, BUFFERSIZE)}
}

// Connect is the method to connect the MSAconverter to some data source
func (proc *GraphSketcher) Connect(previous *MSAconverter) {
	proc.input = previous.output
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *GraphSketcher) Run() {
	defer close(proc.output)

	// after sketching all the received graphs, add the graphs to a store and save it
	graphChan := make(chan *graph.GrootGraph)
	graphStore := make(graph.Store)

	// receive the graphs to be sketched
	var wg sync.WaitGroup
	for newGraph := range proc.input {
		wg.Add(1)
		go func(grootGraph *graph.GrootGraph) {

			if !grootGraph.Masked {

				// create sketch for each window in the graph (merging consecutive windows with identical sketches)
				windows, err := grootGraph.WindowGraph(proc.info.WindowSize, proc.info.KmerSize, proc.info.SketchSize)
				misc.ErrorCheck(err)

				// send the windows on the indexing
				proc.output <- windows
			}

			// this graph is sketched, now send it on to be saved in the current process
			graphChan <- grootGraph
			wg.Done()
		}(newGraph)
	}
	go func() {
		wg.Wait()
		close(graphChan)
	}()

	// collect the graphs
	numMasked := 0
	numWindows := 0
	propDistinctSketches := 0.0
	for sketchedGraph := range graphChan {
		if sketchedGraph.Masked {
			numMasked++
		} else {

			// get the number of windows sketched, the proportion which resulted in distinct sketches, and the max span between merged sketches
			nw, nds, ms, err := sketchedGraph.GetSketchStats()
			misc.ErrorCheck(err)
			propDistinct := float64(nds) / float64(nw)

			// check for the max span between identical sketches
			if ms > proc.info.MaxSketchSpan {
				refs, err := sketchedGraph.GetRefIDs()
				misc.ErrorCheck(err)
				misc.ErrorCheck(fmt.Errorf("graph (ID: %d) encountered where %d sketches in a row were merged (max permitted span: %d)\nencoded seqs: %v", sketchedGraph.GraphID, ms, proc.info.MaxSketchSpan, refs))
			}

			numWindows += nw
			propDistinctSketches += propDistinct
		}

		// store the graph
		graphStore[sketchedGraph.GraphID] = sketchedGraph
	}

	// check some graphs have been sketched
	numGraphs := len(graphStore) - numMasked
	if numGraphs == 0 {
		misc.ErrorCheck(fmt.Errorf("could not create and sketch any graphs"))
	}
	log.Printf("\tnumber of groot graphs built: %d", len(graphStore))
	log.Printf("\t\tgraphs sketched: %d", numGraphs)
	log.Printf("\t\tgraph windows processed: %d", numWindows)
	log.Printf("\t\tmean approximate distinct sketches per graph: %.2f%%", (propDistinctSketches/float64(numGraphs))*100)

	// add the graphs to the pipeline info
	proc.info.Store = graphStore
}

// SketchIndexer is a pipeline process that adds sketches to the LSH Ensemble
type SketchIndexer struct {
	info  *Info
	input chan map[string][]lshe.Key
}

// NewSketchIndexer is the constructor
func NewSketchIndexer(info *Info) *SketchIndexer {
	return &SketchIndexer{info: info}
}

// Connect is the method to connect the MSAconverter to some data source
func (proc *SketchIndexer) Connect(previous *GraphSketcher) {
	proc.input = previous.output
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *SketchIndexer) Run() {

	// create the containment index struct
	numKmers := ((proc.info.WindowSize - proc.info.KmerSize) + 1)
	index := lshe.InitIndex(proc.info.NumPart, proc.info.MaxK, numKmers, proc.info.SketchSize)

	// collect the window sketches from each graph
	sketchCount := 0
	for graphWindowMap := range proc.input {

		// check the windows for each node of the graph
		for keyBase, windows := range graphWindowMap {
			for i, window := range windows {

				// use the iterator to distinguish windows from the same start node
				lookup := fmt.Sprintf("%v-%d", keyBase, i)

				// add to the index
				misc.ErrorCheck(index.AddWindow(lookup, window))
				sketchCount++
			}
		}
	}

	// the index has all the windows, now add it to the runtime info for serialisation
	proc.info.AttachDB(index)
	log.Printf("\tnumber of sketches added to the LSH Ensemble index: %d\n", sketchCount)
}
