package graph

import (
	"fmt"

	"github.com/will-rowe/groot/src/em"
)

// grootGraphPath
type grootGraphPath struct {
	pathID    uint32
	name      []byte
	nodes     []uint64
	sequences [][]byte
	abundance float64
}

// grootGraphPaths is a slice that implements the sort interface, allowing us to sort the paths by abundance
type grootGraphPaths []*grootGraphPath

func (grootGraphPaths grootGraphPaths) Len() int {
	return len(grootGraphPaths)
}
func (grootGraphPaths grootGraphPaths) Swap(i, j int) {
	grootGraphPaths[i], grootGraphPaths[j] = grootGraphPaths[j], grootGraphPaths[i]
}
func (grootGraphPaths grootGraphPaths) Less(i, j int) bool {
	return grootGraphPaths[i].abundance > grootGraphPaths[j].abundance
}

// RunEM is a method to run EM on an approximately weighted variation graph
func (GrootGraph *GrootGraph) RunEM(minIterations, numIterations int) error {

	// make the map of equivalence classes
	ecMap := make(map[uint64][]uint32)
	counts := make(map[uint64]float64)
	for _, node := range GrootGraph.SortedNodes {

		// some nodes will be marked for skipping after pruning, ignore these
		if node.Marked {
			continue
		}
		if _, ok := ecMap[node.SegmentID]; ok {
			return fmt.Errorf("duplicate node ID found in graph")
		}
		ecMap[node.SegmentID] = node.PathIDs
		counts[node.SegmentID] = (float64(node.KmerFreq) / float64(len(node.Sequence)))
	}

	// set up the EM
	em, err := em.NewEM(numIterations, minIterations, GrootGraph.Paths, GrootGraph.Lengths, ecMap, counts)
	if err != nil {
		return err
	}

	// run the EM
	err = em.Run()
	if err != nil {
		return err
	}

	// collect the results
	emIterations, alpha, err := em.Return()

	// update the graph with the results
	GrootGraph.EMiterations = emIterations
	GrootGraph.alpha = alpha
	return err
}

// ProcessEMpaths is a method to process the paths after EM has been run
func (GrootGraph *GrootGraph) ProcessEMpaths(cutoff float64, totalKmers int) error {
	if GrootGraph.EMiterations == 0 {
		return fmt.Errorf("EM has not been run for this graph")
	}
	pathsTotal := 0.0
	rho := make([]float64, len(GrootGraph.alpha))
	for i := 0; i < len(GrootGraph.alpha); i++ {
		pathsTotal += GrootGraph.alpha[i]
	}
	for i := 0; i < len(GrootGraph.alpha); i++ {
		rho[i] = GrootGraph.alpha[i] / pathsTotal
	}

	// rank EM paths, apply cutoff and store
	GrootGraph.abundances = make(map[uint32]float64)
	for i, val := range rho {
		kmerShare := (val * float64(GrootGraph.KmerTotal)) / float64(totalKmers)

		// ensure kept paths are above cutoff
		if kmerShare >= cutoff {
			GrootGraph.abundances[uint32(i)] = kmerShare
			continue
		}
		delete(GrootGraph.Paths, uint32(i))
	}
	return nil
}

// GetEMpaths is a method to print the paths and abundance values
func (GrootGraph *GrootGraph) GetEMpaths() ([]string, []float64) {
	paths := make([]string, len(GrootGraph.grootPaths))
	vals := make([]float64, len(GrootGraph.grootPaths))
	for i, path := range GrootGraph.grootPaths {
		paths[i] = string(path.name)
		vals[i] = path.abundance
	}
	return paths, vals
}

/*
// segmentInfo is a struct used by MCMC to calculate acceptance probabilties
type segmentInfo struct {
	id uint64  // the segmentID field of a node
	tp float64 // the transition probability of arriving at this node
	kf float64 // the k-mer frequency of a node
}

// ByTP implements sort.Interface for []segmentInfo based on the transition probability field
type ByTP []segmentInfo

func (a ByTP) Len() int           { return len(a) }
func (a ByTP) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByTP) Less(i, j int) bool { return a[i].tp > a[j].tp }

// BuildMarkovChain is a method to add the paths from a graph to a markov chain model
func (GrootGraph *GrootGraph) BuildMarkovChain(chain *markov.Chain) error {
	// see if path finding has already been run on this graph
	if len(GrootGraph.mcmcPaths) != 0 {
		return fmt.Errorf("MCMC path finding has already been done on this graph - possible duplicate graph")
	}

	// get all the reference paths
	if len(GrootGraph.originalPaths) == 0 {
		if err := GrootGraph.GetPaths(); err != nil {
			return err
		}
	}

	// add the paths to the markov chain
	for _, path := range GrootGraph.originalPaths {
		// convert sequences from []byte to []string
		seqs := make([]string, len(path.sequences))
		for i, seq := range path.sequences {
			seqs[i] = string(seq)
		}
		chain.Add(seqs)
	}

	// TODO: add some checks, possibly a mark to note that a graph has been added to the chain
	return nil
}

// FindMarkovPaths is a method to collect probable paths through the graph
func (GrootGraph *GrootGraph) FindMarkovPaths(chain *markov.Chain, bootstraps int, scaling float64) error {
	// get the total number of k-mers in this graph
	kmerTotal := 0.0
	for _, node := range GrootGraph.SortedNodes {
		kmerTotal += node.KmerFreq
	}

	// get a copy of the graph's startingNodes
	startingNodes, err := GrootGraph.GetStartNodes()
	if err != nil {
		return err
	}

	// set up to run each bootstrap of the path finder concurrently
	pathSend := make(chan *grootGraphPath)
	errChannel := make(chan error, 1)
	var wg sync.WaitGroup
	go func() {
		wg.Wait()
		close(pathSend)
		close(errChannel)
	}()

	// run all the bootstraps for each starting node
	for i := 0; i < len(startingNodes); i++ {

		// get the starting node, path tracker, memory
		n, pt, pm, err := GrootGraph.setupPathSearch(startingNodes[i], chain.Order)
		if err != nil {
			return err
		}
		wg.Add(bootstraps)

		// run the bootstraps
		for bootstrapCounter := 0; bootstrapCounter < bootstraps; bootstrapCounter++ {

			go func(node *GrootGraphNode, pathTracker grootGraphPath, pathMemory []string) {
				defer wg.Done()

				// start following a path
				for {
					// if there are no more edges in the current path, report the current path and end this bootstrap
					if len(node.OutEdges) == 0 {

						// send the path by a channel
						pathSend <- &pathTracker
						break
					}

					// add the current node sequence to the end of the pathMemory, and remove the oldest item
					pathMemory = append(pathMemory, string(node.Sequence))
					pathMemory = pathMemory[1:]

					// if the current node has only one edge, get ready to follow it
					if len(node.OutEdges) == 1 {
						node, err = GrootGraph.GetNode(node.OutEdges[0])
						if err != nil {
							errChannel <- err
						}
					} else {

						// otherwise, we need to choose which edge we are going to follow
						segmentCandidates := make([]segmentInfo, len(node.OutEdges))

						// start by getting all the weights
						combinedWeights := 0.0
						for i, edgeID := range node.OutEdges {
							edgeNode, err := GrootGraph.GetNode(edgeID)
							if err != nil {
								errChannel <- err
							}
							tp, err := chain.TransitionProbability(string(edgeNode.Sequence), pathMemory)
							if err != nil {
								errChannel <- err
							}
							segmentCandidates[i] = segmentInfo{edgeID, tp, edgeNode.KmerFreq}
							combinedWeights += edgeNode.KmerFreq
						}

						// sort the edges by decreasing transition probability
						sort.Sort(ByTP(segmentCandidates))

						// alternative: should the edge selection be random?
						//rand.Shuffle(len(segmentCandidates), func(i, j int) {
						//	segmentCandidates[i], segmentCandidates[j] = segmentCandidates[j], segmentCandidates[i]
						//})

						// MCMC - metropolis hasting (https://towardsdatascience.com/markov-chain-monte-carlo-291d8a5975ae)
						// generate a random number from a uniform distribution, in the range 0..1
						randomNumber := GrootGraph.rng.Float64Range(0.0, 1.0)
						selected := false
						var selectedEdgeID uint64

						for _, edge := range segmentCandidates {

							// calculate the acceptance probability for this edge
							edgeFreq := (edge.kf / combinedWeights)
							//edgeFreq := 1 - (edge.kf / kmerTotal)
							acceptanceProbability := edgeFreq * edge.tp

							///fmt.Println(edge.kf, combinedWeights, kmerTotal, edgeFreq, edge.tp, acceptanceProbability)

							if acceptanceProbability == math.Inf(1) {
								continue
							}

							// evaluate the acceptance probability against the rng
							if acceptanceProbability > 1.0 || acceptanceProbability >= randomNumber {
								selectedEdgeID = edge.id
								selected = true
								pathTracker.weights = append(pathTracker.weights, acceptanceProbability)
								break
							}

						}

						// check we selected an edge, if we didn't then end the current path search
						if selected == false {
							break
							//errChannel <- fmt.Errorf("could not find suitable edge from node %d", node.SegmentID)
						}

						// get ready to follow the selected node
						node, err = GrootGraph.GetNode(selectedEdgeID)
						if err != nil {
							errChannel <- err
						}
					}

					//if err := chain.Scale(string(node.Sequence), pathMemory, scaling); err != nil {
					//	errChannel <- err
					//}

					// add the ID of the selected node to the tracker before going on to check its edges
					pathTracker.nodes = append(pathTracker.nodes, node.SegmentID)
				}
				return
			}(n, *pt, pm)
		}
	}

	// collect all the paths
	paths := []*grootGraphPath{}
	for path := range pathSend {

		// keep the first path received
		if len(paths) == 0 {
			paths = append(paths, path)
			continue
		}

		// otherwise, check if the path has been seen before
		keep := true
		for _, x := range paths {
			if exists := misc.Uint64SliceEqual(path.nodes, x.nodes); exists {
				keep = false
				break
			}
		}
		if keep {
			pathHolder := ""
			for _, node := range path.nodes {
				pathHolder += fmt.Sprintf("%d+,", node)
			}
			paths = append(paths, path)
		}
	}

	// once all the paths are collected, store them and exit
	GrootGraph.mcmcPaths = paths

	// check there were no errors
	if err := <-errChannel; err != nil {
		return err
	}

	return nil
}

// ProcessMarkovPaths is a method to
func (GrootGraph *GrootGraph) ProcessMarkovPaths(probCutoff float64) error {
	if len(GrootGraph.mcmcPaths) == 0 {
		return fmt.Errorf("no markov paths to process for this graph")
	}
	keptPaths := []*grootGraphPath{}
	for pathCount, path := range GrootGraph.mcmcPaths {
		// get the combined probability of each graph transition in this path
		combProb := path.weights[0]
		for _, prob := range path.weights {
			combProb *= prob
		}
		// discard any paths below the probability threshold
		if combProb < probCutoff {
			continue
		}
		// check this mcmcPath against the original paths
		for _, originalPath := range GrootGraph.originalPaths {
			// use the name if we have it
			if path.pathMatch(originalPath) {
				path.name = originalPath.name
				break
			}
		}
		// if we don't know the name, give it one using the pathCount
		if path.name == nil {
			name := fmt.Sprintf("groot-graph-%d-unknownPath-%d", GrootGraph.GraphID, pathCount)
			path.name = []byte(name)
		}
		// keep this path
		keptPaths = append(keptPaths, path)
	}
	// replace the original paths with only the ones passing the probability theshold
	GrootGraph.mcmcPaths = keptPaths
	if len(GrootGraph.mcmcPaths) == 0 {
		return fmt.Errorf("no markov paths were kept for this graph")
	}
	return nil
}

// GetMarkovPaths is a method to return the identified paths
func (GrootGraph *GrootGraph) GetMarkovPaths() ([]string, error) {
	if len(GrootGraph.mcmcPaths) == 0 {
		return nil, fmt.Errorf("no markov paths were found for this graph")
	}
	// collect each path as a string, in GFA path format
	paths := make([]string, len(GrootGraph.mcmcPaths))
	for i, path := range GrootGraph.mcmcPaths {
		// path holder
		pathHolder := ""
		for _, node := range path.nodes {
			pathHolder += fmt.Sprintf("%d+,", node)
		}
		prefix := fmt.Sprintf("P\t%v\t", string(path.name))
		paths[i] = prefix + pathHolder
	}
	return paths, nil
}

// setupPathSearch is a method used by the path builder that takes a start node id and sets up and return a start node, an empty path tracker and path memory, and any error
func (GrootGraph *GrootGraph) setupPathSearch(nodeID uint64, chainOrder int) (*GrootGraphNode, *grootGraphPath, []string, error) {
	// set up the path tracker
	pathTracker := &grootGraphPath{}

	// add start node ID to the tracker
	pathTracker.nodes = []uint64{nodeID}

	// get the node
	node, err := GrootGraph.GetNode(nodeID)
	if err != nil {
		return nil, nil, nil, err
	}

	// check to make sure this starting node has some edges...
	if len(node.OutEdges) == 0 {
		return nil, nil, nil, fmt.Errorf("graph has an unconnected starting node")
	}

	// create a short term path memory
	pathMemory := make([]string, chainOrder, chainOrder+1)
	for i := 0; i < chainOrder; i++ {
		pathMemory[i] = markov.StartToken
	}
	return node, pathTracker, pathMemory, nil
}

*/
