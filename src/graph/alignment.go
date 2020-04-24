package graph

import (
	"fmt"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/lshforest"
	"github.com/will-rowe/groot/src/seqio"
)

// AlignRead is a method to run a read to graph hierarchical alignment
func (GrootGraph *GrootGraph) AlignRead(read *seqio.FASTQread, mapping *lshforest.Key, references []*sam.Reference) ([]*sam.Record, error) {

	// TODO: move these hardcoded values to CLI options
	MaxClip := 2
	MaxShuffles := 15

	// store the ID of the first node in the seed
	seedNodeID := mapping.Node

	// get the node location in the sorted graph using the lookup map
	nodeLookup, ok := GrootGraph.NodeLookup[seedNodeID]
	if !ok {
		return nil, fmt.Errorf("could not perform node lookup during alignment - possible incorrect seed")
	}

	// run the hierarchical alignment
	IDs := []int{}
	startPos := make(map[int]int)
	startClippedBases := 0
	endClippedBases := 0

	// 1. exact alignment and seed offset shuffling
	var shuffles int
	for shuffles = 0; shuffles <= MaxShuffles; shuffles++ {
		IDs, startPos = GrootGraph.performAlignment(nodeLookup, &read.Seq, int(mapping.OffSet))
		if len(IDs) > 0 {
			break
		}
		mapping.OffSet++
	}

	/*
		TODO:
		// 2. reverse seed shuffling (start at last base of the input edge)
	*/

	// 3. hard clipping the start of the read
	if len(IDs) == 0 {

		// reset the offset
		mapping.OffSet -= uint32(shuffles)

		// make a copy of the sequence for clipping
		clippedSeq := read.Seq
		for i := 1; i <= MaxClip; i++ {
			clippedSeq = clippedSeq[i:]
			IDs, startPos = GrootGraph.performAlignment(nodeLookup, &clippedSeq, int(mapping.OffSet))
			startClippedBases++
			if len(IDs) != 0 {
				break
			}
		}
	}

	// 4. hard clipping the end of the read
	if len(IDs) == 0 {

		// reset the start clip
		startClippedBases = 0

		// make a copy of the sequence for clipping
		clippedSeq := read.Seq
		for i := MaxClip; i > 0; i-- {
			clippedSeq = clippedSeq[:len(clippedSeq)-1]
			IDs, startPos = GrootGraph.performAlignment(nodeLookup, &clippedSeq, int(mapping.OffSet))
			endClippedBases++
			if len(IDs) != 0 {
				break
			}
		}
	}

	// TODO: add alignment step to handle terminal alignments

	// return if no alignments found for this read against this graph
	if len(IDs) == 0 {
		return nil, nil
	}

	// report any alignments
	alignments := []*sam.Record{}
	for alignmentCounter, ID := range IDs {

		// set up the alignment record
		seqLength := len(read.Seq) - endClippedBases - startClippedBases
		record := &sam.Record{
			Name: string(read.ID[1:]),
			Seq:  sam.NewSeq(read.Seq[0:seqLength]),
			Qual: read.Qual[0:seqLength],
		}

		// add the reference
		record.Ref = references[ID]

		// add in the start position for the alignment
		record.Pos = startPos[ID]

		// add a perfect alignment CIGAR for length of exact alignment, plus any hard clipping
		cigar := sam.Cigar{}
		if startClippedBases != 0 {
			cigar = append(cigar, sam.NewCigarOp(sam.CigarHardClipped, startClippedBases))
		}
		cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, seqLength))
		if endClippedBases != 0 {
			cigar = append(cigar, sam.NewCigarOp(sam.CigarHardClipped, endClippedBases))
		}
		record.Cigar = cigar

		// set the MAPQ
		// TODO: this is just left in to have a valid SAM file, I need to set these values correctly
		record.MapQ = 30

		// specify the read orientation and if secondary alignment
		// note: secondary alignments here are just classed as any additional alignment after the first reported
		if len(IDs) > 1 && alignmentCounter != 0 {
			record.Flags |= sam.Secondary
		}
		if read.RC == true {
			record.Flags |= sam.Reverse
		}

		// store the alignment
		alignments = append(alignments, record)
	}

	return alignments, nil
}

// performAlignment does the actual work
func (GrootGraph *GrootGraph) performAlignment(NodeLookup int, read *[]byte, offset int) ([]int, map[int]int) {

	// create some empty variables to store the ID and start Pos of any alignment
	IDs := []int{}
	startPos := make(map[int]int)

	// variables to send and hold paths
	var wg sync.WaitGroup
	sendPath := make(chan []uint64)
	paths := [][]uint64{}
	readLength := len(*read)
	wg.Add(1)
	go func() {
		defer wg.Done()
		_ = GrootGraph.dfsRecursive(GrootGraph.SortedNodes[NodeLookup], read, 0, []uint64{}, sendPath, readLength, offset)
	}()
	go func() {
		wg.Wait()
		close(sendPath)
	}()

	// collect any succecssful paths from the local alignment
	for path := range sendPath {
		paths = append(paths, path)
	}

	// process the traversals
	if len(paths) != 0 {
		IDs, startPos = GrootGraph.processTraversal(paths, offset)
	}
	return IDs, startPos
}

// dfsRecursive is a function to perform an alignment using recursive depth first search of a variation graph
func (GrootGraph *GrootGraph) dfsRecursive(node *GrootGraphNode, read *[]byte, distance int, path []uint64, sendPath chan []uint64, readLength, offset int) bool {

	// check that the offset does not exceed the node sequence length
	if offset >= len(node.Sequence) {
		return false
	}

	// iterate over the segment sequence held by this node and check matches
	for _, base := range node.Sequence[offset:] {

		// stop matching if the read length has been reached
		if distance == readLength {
			break
		}

		// skip to the next base and index position if the reference base is an N. TODO: better handling of non-ACTG bases in reference graphs
		if base == 'N' {
			distance++
			continue
		}

		// increment the distance counter for each match
		if base == (*read)[distance] {
			distance++
		} else {
			return false // terminate this DFS
		}
	}

	// increment the path to include the segment that has just been matched
	path = append(path, node.SegmentID)

	// if we have a consensus length that equals read length (==exact match), or there are no more nodes in the graph - end the DFS and report the alignment path
	if distance == readLength || len(node.OutEdges) == 0 {
		pathCopy := make([]uint64, len(path))
		for i, j := range path {
			pathCopy[i] = j
		}
		sendPath <- pathCopy
		return true
	}

	// a switch to halt the current path traversal upon exact alignment covering whole read
	aligned := false

	// continue along the current traversal by trying to build the consensus using neighbouring nodes
	for _, neighbourNode := range node.OutEdges {
		// get the node from the sorted graph using the lookup
		NodeLookup, ok := GrootGraph.NodeLookup[neighbourNode]
		if !ok {
			panic("could not perform node lookup during alignment - possible incorrect seed")
		}
		// call the DFS func again
		if result := GrootGraph.dfsRecursive(GrootGraph.SortedNodes[NodeLookup], read, distance, path, sendPath, readLength, 0); result == true {
			aligned = true
		}
	}
	return aligned
}

/*
  A function to report a graph traversal relative to a linear reference sequence

  it receives a path(s) (the nodes which were traversed during a successful alignment) and evaluates the parent sequences for each node
  it returns the ID of the most frequently occurring parent sequence, plus the start position of the alignment, relative to the linear reference sequence of the parent(s)
  it also updates the grootGraph to record the aligned read
*/
func (GrootGraph *GrootGraph) processTraversal(paths [][]uint64, offset int) ([]int, map[int]int) {
	IDassignments := []int{}
	startPositions := make(map[int]int)

	// process one or more paths and combine the results
	for pathIterator := 0; pathIterator < len(paths); pathIterator++ {
		nodeIDs := make(map[int]int)
		IDassignment := []int{}
		startPos := make(map[int]int)
		pathLength := len(paths[pathIterator])

		// move along the alignment path
		for i := 0; i < pathLength; i++ {

			// lookup the position of the nodeID in the graph and get the full node info
			NodeLookup, ok := GrootGraph.NodeLookup[paths[pathIterator][i]]
			if !ok {
				panic("could not perform node lookup during alignment - possible incorrect seed")
			}

			// record the alignment against this node NOTE: this has been moved to the weighting code (IncrementSubPath)
			// GrootGraph.SortedNodes[NodeLookup].IncrementReadCount()

			// tally the parent IDs of this node
			node := GrootGraph.SortedNodes[NodeLookup]
			for _, id := range node.PathIDs {
				if _, ok := nodeIDs[int(id)]; ok {
					nodeIDs[int(id)]++
				} else {
					nodeIDs[int(id)] = 1
				}
				// if this is the starting node of the alignment, grab the start position (relative to reference sequence(s))
				if i == 0 {
					startPos[int(id)] = node.Position[int(id)] + offset
				}
			}
		}
		// assign reference to traversal (only reference IDs that are present in each node of the traversal are used)
		for key, value := range nodeIDs {
			if value < pathLength {
				continue
			} else {
				IDassignment = append(IDassignment, key)
			}
		}
		// add the ID for the current path to the return variables
		IDassignments = append(IDassignments, IDassignment...)
		for key, value := range startPos {
			if _, ok := startPositions[key]; !ok {
				startPositions[key] = value
			}
		}
	}
	return IDassignments, startPositions
}
