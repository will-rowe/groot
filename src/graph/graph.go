// Package graph is used to process graphs. It converts, writes, aligns reads and processes GROOT graphs.
package graph

import (
	"bytes"
	"fmt"
	"sort"
	"strconv"
	"sync"

	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/lshforest"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
)

// GrootGraph is the variation graph implementation used by GROOT
type GrootGraph struct {
	GrootVersion        string
	GraphID             uint32
	SortedNodes         []*GrootGraphNode  // essentially, this is the graph - a topologically sorted array of nodes
	Paths               map[uint32][]byte  // lookup to relate PathIDs in each node to a path name
	Lengths             map[uint32]int     // lengths of sequences held in graph (lookup key corresponds to key in Paths)
	NodeLookup          map[uint64]int     // this map returns a the position of a node in the SortedNodes array, using the node segmentID as the locator
	Masked              bool               // a flag to prevent the graph being used by GROOT
	KmerTotal           uint64             // the total number of k-mers projected onto the graph
	EMiterations        int                // the number of EM iterations ran
	alpha               []float64          // indices match the Paths
	abundances          map[uint32]float64 // abundances of kept paths, relative to total k-mers processed during sketching
	grootPaths          grootGraphPaths    // an explicit path through the graph
	numWindows          int                // number of windows that were sketched
	numDistinctSketches int                // number of distinct sketches produced
	maxSpan             uint32             // max span between sketches that have been merged
}

// CreateGrootGraph is a GrootGraph constructor that takes a GFA instance and stores the info as a graph and then runs a topological sort
func CreateGrootGraph(gfaInstance *gfa.GFA, id int) (*GrootGraph, error) {

	// construct an empty graph
	newGraph := &GrootGraph{
		GraphID:             uint32(id),
		Paths:               make(map[uint32][]byte),
		Lengths:             make(map[uint32]int),
		NodeLookup:          make(map[uint64]int),
		numWindows:          0,
		numDistinctSketches: 0,
		maxSpan:             0,
	}

	// collect all the segments from the GFA instance and create the nodes
	segments, err := gfaInstance.GetSegments()
	if err != nil {
		return nil, err
	}
	for nodeIterator, segment := range segments {
		// check the segment name can be stored as an int

		// TODO: will fix the handling of segmentIDs between GFA and GROOT -- need to use uint64
		segID, err := strconv.Atoi(string(segment.Name))
		if err != nil {
			return nil, fmt.Errorf("could not convert segment name from GFA into an int for groot graph: %v", segment.Name)
		}
		// convert all bases to upperCase and check for non-ACTGN chars
		seq := seqio.Sequence{Seq: segment.Sequence}
		if err := seq.BaseCheck(); err != nil {
			return nil, err
		}
		// check if there are optional fields included
		kmerCount := 0.0
		kc, err := segment.GetKmerCount()
		if err != nil {
			return nil, err
		}
		if kc != 0.0 {
			kmerCount = float64(kc)
		}
		newNode := &GrootGraphNode{
			SegmentID:     uint64(segID),
			SegmentLength: float64(len(seq.Seq)),
			Sequence:      seq.Seq,
			KmerFreq:      kmerCount,
			Position:      make(map[int]int),
		}

		// store the new node in the graph and record it's location in the silce by using the NodeLookup map
		newGraph.SortedNodes = append(newGraph.SortedNodes, newNode)
		newGraph.NodeLookup[uint64(segID)] = nodeIterator
		newGraph.KmerTotal += uint64(kmerCount)
	}
	// collect all the links from the GFA instance and add edges to the nodes
	links, err := gfaInstance.GetLinks()
	if err != nil {
		return nil, err
	}
	for _, link := range links {
		// get the from and to segment IDs
		fromSegID, err := strconv.Atoi(string(link.From))
		if err != nil {
			return nil, fmt.Errorf("could not convert fromSegID name from GFA into an int for groot graph: %v", link.From)
		}
		toSegID, err := strconv.Atoi(string(link.To))
		if err != nil {
			return nil, fmt.Errorf("could not convert toSegID name from GFA into an int for groot graph: %v", link.To)
		}
		// add the outEdges
		nodeLocator := newGraph.NodeLookup[uint64(fromSegID)]
		newGraph.SortedNodes[nodeLocator].OutEdges = append(newGraph.SortedNodes[nodeLocator].OutEdges, uint64(toSegID))
	}
	// collect all the paths from the GFA instance and add pathIDs to each node
	// this function also populates the reference Position for each node
	paths, err := gfaInstance.GetPaths()
	if err != nil {
		return nil, err
	}
	for pathIterator := uint32(0); pathIterator < uint32(len(paths)); pathIterator++ {
		// add the path name to the lookup
		newGraph.Paths[pathIterator] = paths[pathIterator].PathName
		for _, seg := range paths[pathIterator].SegNames {
			// strip the plus
			seg = bytes.TrimSuffix(seg, []byte("+"))
			segID, err := strconv.Atoi(string(seg))
			if err != nil {
				return nil, fmt.Errorf("could not convert segment name from GFA path into an int for groot graph: %v\n%v", string(seg), string(paths[pathIterator].PathName))
			}
			nodeLocator := newGraph.NodeLookup[uint64(segID)]
			newGraph.SortedNodes[nodeLocator].PathIDs = append(newGraph.SortedNodes[nodeLocator].PathIDs, pathIterator)
			// add the first segment of this path to the start nodes
			//if i == 0 {
			//	newGraph.startNodes[uint64(segID)] = struct{}{}
			//}
		}
	}
	// return without toposort if only one node present (graph with single sequence)
	if len(newGraph.SortedNodes) > 1 {
		err = newGraph.topoSort()
	}
	// get and store the lengths of each sequence held in the graph
	seqs, err := newGraph.Graph2Seqs()
	if err != nil {
		return nil, err
	}
	for pathID, path := range seqs {
		newGraph.Lengths[uint32(pathID)] = len(path)
	}
	// return the new GrootGraph
	return newGraph, err
}

// topoSort runs a topological sort on the GrootGraph
func (GrootGraph *GrootGraph) topoSort() error {
	// copy all of the graph nodes into a map (so we can keep track of what we have processed)
	nodeMap := make(map[uint64]*GrootGraphNode)
	toposortStart := []uint64{}
	seenPaths := make(map[uint32]int)
	for _, node := range GrootGraph.SortedNodes {
		if len(seenPaths) == len(GrootGraph.Paths) {
			break
		}
		// if this node is from a sequence we have not seen yet, mark the node as a starting node for toposort
		for _, path := range node.PathIDs {
			if _, ok := seenPaths[path]; !ok {
				toposortStart = append(toposortStart, node.SegmentID)
			}
		}
		// check for duplicate nodes
		if _, ok := nodeMap[node.SegmentID]; ok {
			return fmt.Errorf("graph contains duplicate nodes (identical segment IDs)")
		}
		// add the node to a map
		nodeMap[node.SegmentID] = node
	}
	// clear the SortedNodes and NodeLookup
	GrootGraph.SortedNodes = []*GrootGraphNode{}
	GrootGraph.NodeLookup = make(map[uint64]int)
	// run the topological sort  - try starting from each node that was in the first slot of the nodeholder (start of the MSA)
	seen := make(map[uint64]struct{})
	for len(nodeMap) > 1 {
		for _, start := range toposortStart {
			if _, ok := nodeMap[start]; !ok {
				continue
			}
			GrootGraph.traverse(nodeMap[start], nodeMap, seen)
		}
	}
	// check all traversals have been taken
	if len(nodeMap) > 0 {
		return fmt.Errorf("topological sort failed - too many nodes remaining in the pre-sort list")
	}
	return nil
}

//  traverse is a helper method to topologically sort a graph
func (GrootGraph *GrootGraph) traverse(node *GrootGraphNode, nodeMap map[uint64]*GrootGraphNode, seen map[uint64]struct{}) {
	// skip if we are already handling the current node
	if _, ok := seen[node.SegmentID]; ok {
		return
	}
	// make sure node is still in the graph
	if _, ok := nodeMap[node.SegmentID]; ok {
		// record that we are processing this node
		seen[node.SegmentID] = struct{}{}
		// sort the output nodes for this node and then traverse them in reverse order
		sort.Sort(sort.Reverse(node.OutEdges))
		for i, j := 0, len(node.OutEdges); i < j; i++ {
			// check if the outedges have been traversed
			if _, ok := nodeMap[node.OutEdges[i]]; !ok {
				continue
			}
			GrootGraph.traverse(nodeMap[node.OutEdges[i]], nodeMap, seen)
		}
		// delete the node from the temporary holders
		delete(nodeMap, node.SegmentID)
		delete(seen, node.SegmentID)
		// update the sorted node slice and the lookup
		GrootGraph.SortedNodes = append([]*GrootGraphNode{node}, GrootGraph.SortedNodes...)
		GrootGraph.NodeLookup[node.SegmentID] = len(nodeMap)
	}
}

// GetSketchStats returns the number of windows processed from the graph, the number of distinct sketches, the max span for merged sketches and any error
func (GrootGraph *GrootGraph) GetSketchStats() (int, int, int, error) {
	if GrootGraph.numWindows == 0 {
		return 0, 0, 0, fmt.Errorf("graph has not been sketched yet")
	}
	return GrootGraph.numWindows, GrootGraph.numDistinctSketches, int(GrootGraph.maxSpan), nil
}

// WindowGraph is a method to slide a window over each path through the graph, sketching the paths and getting window information
func (GrootGraph *GrootGraph) WindowGraph(windowSize, kmerSize, sketchSize int) (map[string][]lshforest.Key, error) {

	// get the linear sequences for this graph
	pathSeqs, err := GrootGraph.Graph2Seqs()
	if err != nil {
		return nil, err
	}

	// update the graph with the number of windows to expect
	GrootGraph.numWindows = 0
	for _, len := range GrootGraph.Lengths {
		GrootGraph.numWindows += len - windowSize + 1
	}

	// window each path
	pathWindows := make(chan lshforest.Key, 100)
	var pathWG sync.WaitGroup
	pathWG.Add(len(GrootGraph.Paths))
	for pathID := range GrootGraph.Paths {

		// get the length of the linear reference for this path
		pathLength := GrootGraph.Lengths[pathID]

		// add another debug panic here in case being asked to window graph with short seqs
		if pathLength < windowSize {
			return nil, fmt.Errorf("graph contains sequence < window size")
		}

		// get the sequence for this path
		pathSequence := pathSeqs[pathID]

		// process each path in a go routine, use a copy of the nodes to prevent races
		go func(sortedNodes []*GrootGraphNode, pathID uint32, pathSequence []byte, pathLength int) {
			defer pathWG.Done()

			// for each base in the linear reference sequence, get the segmentID and offset of its location in the graph
			segs := make([]uint64, pathLength, pathLength)
			offSets := make([]uint32, pathLength, pathLength)
			iterator := 0
			for _, node := range sortedNodes {
				for _, id := range node.PathIDs {
					if id == pathID {
						for offset := uint32(0); offset < uint32(len(node.Sequence)); offset++ {
							segs[iterator] = node.SegmentID
							offSets[iterator] = offset
							iterator++
						}
					}
				}
			}
			if iterator != pathLength {
				panic("windowing did not traverse entire path")
			}

			// hold a window until a new sketch is encountered
			var windowHolder lshforest.Key
			sketchSent := false

			// start windowing the path sequence
			numWindows := pathLength - windowSize + 1
			for i := 0; i < numWindows; i++ {

				// sketch the current window
				windowSeq := seqio.Sequence{Seq: pathSequence[i : i+windowSize]}
				sketch, err := windowSeq.RunMinHash(kmerSize, sketchSize, false, nil)
				if err != nil {
					panic(err)
				}

				// if this is not the first window, check if the current sketch matches the previous sketch and window
				merge := false
				if i != 0 {

					// if sketch doesn't match previous we send the old window on, otherwise we merge current window into previous one
					if !misc.Uint64SliceEqual(windowHolder.Sketch, sketch) {
						pathWindows <- windowHolder
						sketchSent = true
					} else {
						merge = true
					}
				}

				// if the first window, or we have just sent a window on, init a windowHolder
				if !merge {
					windowHolder = lshforest.Key{
						GraphID:        GrootGraph.GraphID,
						Node:           segs[i],
						OffSet:         offSets[i],
						ContainedNodes: make(map[uint64]float64),
						Ref:            []uint32{pathID},
						Sketch:         sketch,
						MergeSpan:      0,
						WindowSize:     uint32(windowSize),
					}
				}

				// regardless of merge, need to add current windows nodes to the windowHolder's map
				for _, y := range segs[i : i+windowSize] {
					windowHolder.ContainedNodes[uint64(y)]++
				}

				// if merging, update the merge span (number of consecutive windows with same sketch)
				if merge {
					windowHolder.MergeSpan++
				}

				// if we've got to the final window and no sketches have been sent (i.e. long window or looong merge), send the sketch
				if !sketchSent && i == (numWindows-1) {
					pathWindows <- windowHolder
				}
			}
		}(GrootGraph.SortedNodes, pathID, pathSequence, pathLength)
	}

	// wait and close the path window channel
	go func() {
		pathWG.Wait()
		close(pathWindows)
	}()

	// collect sketched windows from all paths and merge identical windows from different paths if same start node+offset
	windowLookup := make(map[string][]lshforest.Key)
	for window := range pathWindows {

		// convert the graph window data to a key that links the sketch to the graphID, start node and offset
		key := fmt.Sprintf("g%dn%do%d", window.GraphID, window.Node, window.OffSet)

		// use the key to check if the there are multiple windows in this graph from the same node+offset
		if existingWindowLocation, ok := windowLookup[key]; ok {

			// check for duplicate sketches at the same node+offset
			duplicateSketch := false
			for _, existingWindow := range existingWindowLocation {

				// if the sketches match, merge the window into the existing one
				if misc.Uint64SliceEqual(existingWindow.GetSketch(), window.Sketch) {
					for node, freq := range window.ContainedNodes {
						existingWindow.ContainedNodes[node] += freq
					}
					existingWindow.Ref = append(existingWindow.Ref, window.Ref...)

					// keep the greatest merge span from any consecutive same-path merges
					if window.MergeSpan > existingWindow.MergeSpan {
						existingWindow.MergeSpan = window.MergeSpan
					}
					duplicateSketch = true
					break
				}
			}

			// if not a duplicate sketch, add the new window to the current window location
			if !duplicateSketch {
				windowLookup[key] = append(windowLookup[key], window)
				GrootGraph.numDistinctSketches++
			}
		} else {
			windowLookup[key] = []lshforest.Key{window}
			GrootGraph.numDistinctSketches++
		}
	}

	// check we've got some sketches
	if GrootGraph.numDistinctSketches == 0 {
		r, _ := GrootGraph.GetRefIDs()
		return nil, fmt.Errorf("no sketches produced after windowing graph seqs: %v", r)
	}
	return windowLookup, nil
}

// IncrementSubPath is a method to adjust the weight of segments that are contained within a given sketch
//  - given a ContainedNodes through a graph, the offset in the first segment, and the number of k-mers in this ContainedNodes
//  - increment the weight of each segment contained in that ContainedNodes by their share of the k-mer coverage for the sketch
func (GrootGraph *GrootGraph) IncrementSubPath(ContainedNodes map[uint64]float64, numKmers float64) error {

	// check the ContainedNodes contains segments
	if len(ContainedNodes) < 1 {
		return fmt.Errorf("ContainedNodes encountered that does not include any segments")
	}

	// if the ContainedNodes is only one segment, then it is straightforward to increment
	if len(ContainedNodes) == 1 {
		for nodeID := range ContainedNodes {
			node, err := GrootGraph.GetNode(nodeID)
			if err != nil {
				return fmt.Errorf("could not perform nodelookup to increment ContainedNodes weight")
			}

			// give this segment all the k-mers for this sketch
			if err := node.IncrementKmerFreq(numKmers); err != nil {
				return err
			}
		}
		return nil
	}

	// otherwise, there are multiple segments in the path and we now work out the proportion of k-mer coverage belonging to each segment

	// single pass to get total length of ContainedNodes
	totalLength := 0.0
	for segment := range ContainedNodes {
		node, err := GrootGraph.GetNode(segment)
		if err != nil {
			return err
		}
		totalLength += node.SegmentLength
	}

	// now iterate over the segments again, incrementing their weights
	for segment, segCount := range ContainedNodes {
		node, err := GrootGraph.GetNode(segment)
		if err != nil {
			return err
		}
		kmerShare := ((node.SegmentLength / totalLength) * numKmers) * float64(segCount)
		if err := node.IncrementKmerFreq(kmerShare); err != nil {
			return err
		}
	}

	// record the number of kmers projected onto the graph
	GrootGraph.IncrementKmerCount(uint64(numKmers))
	return nil
}

// Prune is a method to remove paths and segments from the graph if they have insufficient coverage
// returns false if pruning results in no paths through the graph remaining
func (GrootGraph *GrootGraph) Prune(minKmerCoverage float64) bool {
	removePathID := make(map[uint32]struct{})
	removeNode := make(map[uint64]struct{})

	// first pass through the graph
	for _, node := range GrootGraph.SortedNodes {

		// get the per-base k-mer frequency for the node
		perbaseCoverage := node.KmerFreq / node.SegmentLength

		// if the coverage is too low for this node, add the segmentID and the contained pathIDs to the removal lists
		if perbaseCoverage < minKmerCoverage {
			for _, id := range node.PathIDs {
				removePathID[id] = struct{}{}
				removeNode[node.SegmentID] = struct{}{}
			}
		}
	}

	// if all the paths need removing, just exit now!
	if len(removePathID) == len(GrootGraph.Paths) {
		return false
	}

	// if it doesn't need pruning, return true
	if len(removeNode) == 0 {
		return true
	}

	// second pass through the graph to prune all the marked nodes and paths
	// TODO: I'm just creating a new slice at the moment and copying nodes which aren't marked
	// TODO: shall I try popping elements out of the original slices instead -- is that more efficient?
	for i, node := range GrootGraph.SortedNodes {

		// remove marked paths
		updatedPathIDs := make([]uint32, 0, len(node.PathIDs))
		for _, id := range node.PathIDs {
			if _, marked := removePathID[id]; !marked {
				updatedPathIDs = append(updatedPathIDs, id)
			}
		}
		node.PathIDs = updatedPathIDs

		// delete any marked nodes
		if _, marked := removeNode[node.SegmentID]; marked {

			// TODO: I've set the node to nil in the sorted node array - in order to keep the NodeLookup in order. But this isn't pretty and now requires you to check for nil when using the SortedNodes array
			GrootGraph.SortedNodes[i].Marked = true
			delete(GrootGraph.NodeLookup, node.SegmentID)
		}

		// remove any edges referencing deleted nodes
		updatedEdges := make([]uint64, 0, len(node.OutEdges))
		for _, edge := range node.OutEdges {
			if _, marked := removeNode[edge]; !marked {
				updatedEdges = append(updatedEdges, edge)
			}
		}
		node.OutEdges = updatedEdges
	}

	// if a path was removed by pruning, set it's length to 0
	for id := range removePathID {
		if _, path := GrootGraph.Paths[id]; path {

			//delete(GrootGraph.Paths, id)
			GrootGraph.Lengths[id] = 0
		}
	}
	return true
}

// GetNode takes a nodeID and returns a pointer to the corresponding node struct in the graph
// Note: this is method does not allow concurrent map access
func (GrootGraph *GrootGraph) GetNode(nodeID uint64) (*GrootGraphNode, error) {
	NodeLookup, ok := GrootGraph.NodeLookup[nodeID]
	if !ok {
		return nil, fmt.Errorf("can't find node %d in graph", nodeID)
	}
	return GrootGraph.SortedNodes[NodeLookup], nil
}

/*
// GetStartNodes is a method to return a slice of all the node ids which are the first node in a path
func (GrootGraph *GrootGraph) GetStartNodes() ([]uint64, error) {
	if len(GrootGraph.startNodes) == 0 {
		return nil, fmt.Errorf("this graph has no paths")
	}
	// convert the startingNodes from map keys into a slice
	startingNodes := []uint64{}
	for i := range GrootGraph.startNodes {
		startingNodes = append(startingNodes, i)
	}
	if len(startingNodes) == 0 {
		return nil, fmt.Errorf("this graph has no paths")
	}
	return startingNodes, nil
}
*/

// RemoveDeadPaths is a method to remove pathIDs from nodes if the path is no longer present in the graph
func (GrootGraph *GrootGraph) RemoveDeadPaths() error {
	for _, node := range GrootGraph.SortedNodes {

		// some nodes will be marked for skipping after pruning, ignore these
		if node.Marked {
			continue
		}
		updatedPathIDs := []uint32{}
		for _, pathID := range node.PathIDs {
			if _, ok := GrootGraph.Paths[pathID]; ok {
				updatedPathIDs = append(updatedPathIDs, pathID)
			}
		}
		node.PathIDs = updatedPathIDs
	}
	return GrootGraph.GetPaths()
}

// GetPaths is a method to get the paths from a graph
func (GrootGraph *GrootGraph) GetPaths() error {
	if len(GrootGraph.Paths) == 0 {
		return fmt.Errorf("no paths recorded in current graph")
	}
	if GrootGraph.abundances == nil {
		GrootGraph.abundances = make(map[uint32]float64)
	}
	GrootGraph.grootPaths = make(grootGraphPaths, len(GrootGraph.Paths))
	counter := 0

	// TODO: does it make the most sense to loop this way? I.E. loop entire node set for each path?
	for pathID, pathName := range GrootGraph.Paths {
		segIDs := []uint64{}
		segSeqs := [][]byte{}
		refLength := 0
		for _, node := range GrootGraph.SortedNodes {

			// some nodes will be marked for skipping after pruning, ignore these
			//if node.Marked {
			//	continue
			//}

			// if this node is in the path, get the deets
			for _, id := range node.PathIDs {
				if id == pathID {

					// update the Position field of the current node with the current path
					node.Position[int(pathID)] = refLength
					refLength += len(node.Sequence)

					// build the path
					segIDs = append(segIDs, node.SegmentID)
					segSeqs = append(segSeqs, node.Sequence)
				}
			}
		}
		if _, ok := GrootGraph.abundances[pathID]; !ok {
			GrootGraph.abundances[pathID] = 0.0
		}

		// store this path
		GrootGraph.grootPaths[counter] = &grootGraphPath{pathID: pathID, name: pathName, nodes: segIDs, sequences: segSeqs, abundance: GrootGraph.abundances[pathID]}
		counter++
	}
	// sort the paths
	sort.Sort(GrootGraph.grootPaths)
	return nil
}

// Graph2Seqs is a method to convert a variation graph to linear reference sequences
func (GrootGraph *GrootGraph) Graph2Seqs() (map[uint32][]byte, error) {

	// get the paths
	if err := GrootGraph.GetPaths(); err != nil {
		return nil, err
	}

	// create the map - the keys link the sequence to pathID
	seqs := make(map[uint32][]byte)

	// for each path, combine the segment sequences, add to the map
	for _, path := range GrootGraph.grootPaths {
		newSeq := []byte{}
		for i := 0; i < len(path.sequences); i++ {
			newSeq = append(newSeq, path.sequences[i]...)
		}
		seqs[path.pathID] = newSeq
	}
	return seqs, nil
}

// GetRefIDs is a method to return all the reference IDs that are encoded as sequences in the graph
func (GrootGraph *GrootGraph) GetRefIDs() ([]string, error) {
	if len(GrootGraph.Paths) == 0 {
		return nil, fmt.Errorf("no paths registered in graph")
	}
	paths := make([]string, len(GrootGraph.Paths))
	for i, path := range GrootGraph.grootPaths {
		paths[i] = string(path.name)
	}
	return paths, nil
}

// IncrementKmerCount is a method to increment the counter for the number of kmers projected onto the graph
func (GrootGraph *GrootGraph) IncrementKmerCount(increment uint64) {
	GrootGraph.KmerTotal += increment
}
