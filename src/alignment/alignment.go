package alignment

import (
	"fmt"
	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
	"sync"
)

/*
  A function to perform hierarchical local alignment of a read against a graph
*/
func Align(read seqio.FASTQread, seedID int, graph *graph.GrootGraph, refs []*sam.Reference, maxClip int) <-chan *sam.Record {
	// store the ID of the first node in the seed
	seedNodeID := read.Seeds[seedID].Node
	// get the node location in the sorted graph using the lookup map
	NodeLookup, ok := graph.NodeLookup[seedNodeID]
	if !ok {
		misc.ErrorCheck(fmt.Errorf("could not perform node lookup during alignment - possible incorrect seed"))
	}
	// get the offset in the node
	offset := read.Seeds[seedID].OffSet
	// create the return channel
	returnAlignments := make(chan *sam.Record)
	// run the hierarchical alignment
	go func() {
		defer close(returnAlignments)
		IDs := []int{}
		startPos := make(map[int]int)
		// try exact alignment
		IDs, startPos = performAlignment(NodeLookup, &read.Seq, graph, offset, 0)
		// if unsuccessful then try shuffling the seed forward
		if len(IDs) != 0 {
			shuffles := 5
			for shuffles > 0 {
				//for nodeShift < 5 {
				IDs, startPos = performAlignment(NodeLookup, &read.Seq, graph, offset, shuffles)
				if len(IDs) != 0 {
					break
				}
				shuffles--
			}
		}
		// if still unsuccessful, try clipping the end of the read
		hardClip := 0
		clippedSeq := read.Seq
		if len(IDs) == 0 {
			for i := maxClip; i > 0; i-- {
				clippedSeq = clippedSeq[:len(clippedSeq)-1]
				IDs, startPos = performAlignment(NodeLookup, &clippedSeq, graph, offset, 0)
				hardClip++
				if len(IDs) != 0 {
					break
				}
			}
		}
		// TODO: add alignment step to handle terminal alignments
		// END hierarchical alignment
		// report any alignments
		if len(IDs) != 0 {
			for _, ID := range IDs {
				record := &sam.Record{
					Name: string(read.ID[1:]),
					Seq:  sam.NewSeq(read.Seq[0:len(read.Seq)-hardClip]),
					Qual: read.Qual[0:len(read.Seq)-hardClip],
				}
				// add the reference
				record.Ref = refs[ID]
				// add in the start position for the alignment
				record.Pos = startPos[ID]
				// add a perfect alignment CIGAR for length of exact alignment
				cigar := sam.Cigar{}
				if hardClip != 0 {
					cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, len(clippedSeq)))
					cigar = append(cigar, sam.NewCigarOp(sam.CigarHardClipped, hardClip))
				} else {
					cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, len(read.Seq)))
				}
				record.Cigar = cigar
				record.MapQ = 30 // TODO: this is just left in to have a valid SAM file, I need to set these values correctly
				// TODO: specify the read orientation and if secondary alignment
				// send the record back
				returnAlignments <- record
			}
		}
	}()
	return returnAlignments
}

func performAlignment(NodeLookup int, read *[]byte, graph *graph.GrootGraph, offset, shuffle int) ([]int, map[int]int) {
	// create some empty variables to store the ID and start Pos of any alignment
	IDs := []int{}
	startPos := make(map[int]int)
	// variables to send and hold paths
	var wg sync.WaitGroup
	sendPath := make(chan []int)
	paths := [][]int{}
	readLength := len(*read)
	wg.Add(1)
	go func() {
		defer wg.Done()
		_ = DFSrecursive(graph.SortedNodes[NodeLookup], graph, read, 0, []int{}, sendPath, readLength, offset, shuffle)
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
		IDs, startPos = processTraversal(graph, paths, offset)
	}
	return IDs, startPos
}

/*
  A function to perform an alignment using recursive depth first search of a variation graph
*/
func DFSrecursive(node *graph.GrootGraphNode, graph *graph.GrootGraph, read *[]byte, distance int, path []int, sendPath chan []int, readLength, offset, shuffle int) bool {
	// a switch to halt the current path traversal upon exact alignment covering whole read
	aligned := false
	// iterate over the segment sequence held by this node and check matches
	for _, base := range node.Sequence[offset:] {
		// stop matching if the read length has been reached
		if distance == readLength {
			break
		}
		// increment the distance counter for each match
		if base == (*read)[distance] {
			distance++
		// if no match but performing shuffled alignment, decrement the shuffle
		} else if (len(path)== 0) && (shuffle > 0) {
				shuffle--
				continue
		// if no shuffle or shuffle is finished, terminate this DFS
		} else {
				return false
		}
	}
	// increment the path to include the segment that has just been matched
	path = append(path, node.SegmentID)
	// if we have a consensus length that equals read length (==exact match), or there are no more nodes in the graph - end the DFS and report the alignment path
	if distance == readLength || len(node.OutEdges) == 0 {
		pathCopy := make([]int, len(path))
		for i, j := range path {
			pathCopy[i] = j
		}
		sendPath <- pathCopy
		return true
	}
	// continue along the current traversal by trying to build the consensus using neighbouring nodes
	for _, neighbourNode := range node.OutEdges {
		// get the node from the sorted graph using the lookup
		NodeLookup, ok := graph.NodeLookup[neighbourNode]
		if !ok {
			misc.ErrorCheck(fmt.Errorf("could not perform node lookup during alignment - possible incorrect seed"))
		}
		// call the DFS func again
		if result := DFSrecursive(graph.SortedNodes[NodeLookup], graph, read, distance, path, sendPath, readLength, 0, shuffle); result == true {
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
func processTraversal(graph *graph.GrootGraph, paths [][]int, offset int) ([]int, map[int]int) {
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
			NodeLookup, ok := graph.NodeLookup[paths[pathIterator][i]]
			if !ok {
				misc.ErrorCheck(fmt.Errorf("could not perform node lookup during alignment - possible incorrect seed"))
			}
			node := graph.SortedNodes[NodeLookup]
			// record the alignment against this node
			node.IncrementReadCount()
			// tally the parent IDs of this node
			for _, id := range node.PathIDs {
				if _, ok := nodeIDs[id]; ok {
					nodeIDs[id]++
				} else {
					nodeIDs[id] = 1
				}
				// if this is the starting node of the alignment, grab the start position (relative to reference sequence(s))
				if i == 0 {
					startPos[id] = node.Position[id] + offset
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
