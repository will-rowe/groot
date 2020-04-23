package pipeline

/*
 this part of the pipeline will read in the weighted graphs from the sketch command, find paths in them and return the haplotypes
*/

import (
	"log"
	"regexp"
	"strconv"
	"sync"

	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/version"
)

// GFAreader is a pipeline process that reads in the weighted GFAs
type GFAreader struct {
	info   *Info
	input  []string
	output chan *graph.GrootGraph
}

// NewGFAreader is the constructor
func NewGFAreader(info *Info) *GFAreader {
	return &GFAreader{info: info, output: make(chan *graph.GrootGraph, BUFFERSIZE)}
}

// Connect is the method to connect the GFAreader to some data source
func (proc *GFAreader) Connect(input []string) {
	proc.input = input
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *GFAreader) Run() {
	var wg sync.WaitGroup
	for i, gfaFile := range proc.input {
		gfaObj, err := graph.LoadGFA(gfaFile)
		misc.ErrorCheck(err)

		// grab the total k-mer count across all the graphs
		if i == 0 {
			commentLines := gfaObj.PrintComments()
			re := regexp.MustCompile(`graphs: (\d+)\)`)
			matches := re.FindStringSubmatch(commentLines)
			kmerCount, err := strconv.Atoi((matches[1]))
			misc.ErrorCheck((err))
			proc.info.Haplotype.TotalKmers = kmerCount
		}

		// convert GFAs to GrootGraph and send them on to the path finder
		wg.Add(1)
		go func(gfaID int, g *gfa.GFA) {
			defer wg.Done()
			grootGraph, err := graph.CreateGrootGraph(g, gfaID)
			if err != nil {
				log.Fatal(err)
			}
			proc.output <- grootGraph
		}(i, gfaObj)
	}
	wg.Wait()
	close(proc.output)
}

// EMpathFinder is a pipeline process to identify graph paths using Expectation Maximization
type EMpathFinder struct {
	info   *Info
	input  chan *graph.GrootGraph
	output chan *graph.GrootGraph
}

// NewEMpathFinder is the constructor
func NewEMpathFinder(info *Info) *EMpathFinder {
	return &EMpathFinder{info: info, output: make(chan *graph.GrootGraph)}
}

// Connect is the method to connect the MCMCpathFinder to the output of a GFAreader
func (proc *EMpathFinder) Connect(previous *GFAreader) {
	proc.input = previous.output
}

// ConnectPruner is the method to connect the MCMCpathFinder to the output of a GraphPruner
func (proc *EMpathFinder) ConnectPruner(previous *GraphPruner) {
	proc.input = previous.output
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *EMpathFinder) Run() {
	var wg sync.WaitGroup

	// collect the weighted graphs
	for inputGraph := range proc.input {
		wg.Add(1)

		// add each graph to the graphStore
		proc.info.Store[inputGraph.GraphID] = inputGraph

		// concurrently process the graphs
		go func(g *graph.GrootGraph) {
			defer wg.Done()

			// remove dead ends
			misc.ErrorCheck(g.RemoveDeadPaths())

			// run the EM
			err := g.RunEM(proc.info.Haplotype.MinIterations, proc.info.Haplotype.MaxIterations)
			misc.ErrorCheck(err)

			// process the EM results
			misc.ErrorCheck(g.ProcessEMpaths(proc.info.Haplotype.Cutoff, proc.info.Haplotype.TotalKmers))

			// send the graph to the next process
			proc.output <- g
		}(inputGraph)
	}
	wg.Wait()
	close(proc.output)
}

// HaplotypeParser is a pipeline process to parse the paths produced by the MCMCpathFinder process
type HaplotypeParser struct {
	info   *Info
	input  chan *graph.GrootGraph
	output []string
}

// NewHaplotypeParser is the constructor
func NewHaplotypeParser(info *Info) *HaplotypeParser {
	return &HaplotypeParser{info: info}
}

// Connect is the method to connect the HaplotypeParser to the output of a EMpathFinder
func (proc *HaplotypeParser) Connect(previous *EMpathFinder) {
	proc.input = previous.output
}

// CollectOutput is a method to return what paths are found via MCMC
func (proc *HaplotypeParser) CollectOutput() []string {
	return proc.output
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *HaplotypeParser) Run() {
	meanEMiterations := 0
	keptGraphs := make(graph.Store)
	keptPaths := []string{}
	for g := range proc.input {
		meanEMiterations += g.EMiterations

		// check graph has some paths left
		if len(g.Paths) == 0 {
			continue
		}

		// remove dead ends
		misc.ErrorCheck(g.RemoveDeadPaths())

		// print some stuff
		paths, abundances := g.GetEMpaths()
		log.Printf("\tgraph %d has %d called alleles after EM", g.GraphID, len(paths))
		for i, path := range paths {
			log.Printf("\t- [%v (abundance: %.3f)]", path, abundances[i])
			keptPaths = append(keptPaths, path)
		}
		g.GrootVersion = version.VERSION
		keptGraphs[g.GraphID] = g
	}
	proc.info.Store = keptGraphs
	proc.output = keptPaths
	if len(keptGraphs) == 0 {
		return
	}
	log.Printf("summarising...")
	log.Printf("\tmean number of EM iterations: %d\n", meanEMiterations/len(keptGraphs))
	log.Printf("\tnumber of graphs with viable paths: %d\n", len(keptGraphs))
	log.Printf("\tnumber of called alleles: %d\n", len(keptPaths))

}

/*
// MCMCpathFinder is a pipeline process to identify graph paths using MCMC
type MCMCpathFinder struct {
	info   *Info
	input  chan *graph.GrootGraph
	output chan *graph.GrootGraph
}

// NewMCMCpathFinder is the constructor
func NewMCMCpathFinder(info *Info) *MCMCpathFinder {
	return &MCMCpathFinder{info: info, output: make(chan *graph.GrootGraph)}
}

// Connect is the method to connect the MCMCpathFinder to the output of a GFAreader
func (proc *MCMCpathFinder) Connect(previous *GFAreader) {
	proc.input = previous.output
}

// Run is the method to run this process, which satisfies the pipeline interface
func (proc *MCMCpathFinder) Run() {
	chain := markov.NewChain(proc.info.Haplotype.ChainOrder)
	for graph := range proc.input {
		// add the graph to the graphStore
		proc.info.Store[graph.GraphID] = graph
		// build the markov model
		err := graph.BuildMarkovChain(chain)
		misc.ErrorCheck(err)
	}

	// once all the graphs have been received and added to the combined model, start on the path finding
	var wg sync.WaitGroup
	for _, storedGraph := range proc.info.Store {
		wg.Add(1)
		go func(g *graph.GrootGraph, c *markov.Chain) {
			defer wg.Done()

				err := g.FindMarkovPaths(c, proc.info.Haplotype.BootStraps, proc.info.Haplotype.ScalingFactor)
				misc.ErrorCheck(err)

				err = g.ProcessMarkovPaths(proc.info.Haplotype.ProbabilityThreshold)
				misc.ErrorCheck(err)

				paths, err := g.GetMarkovPaths()
				misc.ErrorCheck(err)

				log.Printf("\tgraph %d has %d markov paths passing thresholds", g.GraphID, len(paths))
				for _, path := range paths {
					log.Printf("\t- [%v]", path)
				}

				// replace the old paths with the new MCMC derived paths
				err = g.PathReplace()
				misc.ErrorCheck(err)

			// send the graph to the next process
			proc.output <- g
		}(storedGraph, chain)
	}
	wg.Wait()
	close(proc.output)
}
*/
