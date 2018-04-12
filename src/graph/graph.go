// this package is used to process GFA formatted graphs and index them for groot
package graph

import (
	"bytes"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"os"
	"sort"
	"strconv"
	"sync"
	"time"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/seqio"
)

/*
   GrootGraphNode is a GFA segment (plus the extra info from path, links etc.)
*/
type GrootGraphNode struct {
	SegmentID int
	Sequence  []byte
	OutEdges  []int
	PathIDs   []int       // a reference to the paths that use this segment (equivalent to the linear reference from which this segment is derived) (value corresponds to key in GrootGraph.Paths)
	Position  map[int]int // the start position of this segment in each reference sequence, using 1-based indexing (lookup key corresponds to key in GrootGraph.Paths)
	ReadCount int         // used during the alignment to record the number of reads covering a node
}

// IncrementReadCount uses a mutex to increase the read count of a node safely during alignment
func (node *GrootGraphNode) IncrementReadCount() {
	var mutex = &sync.Mutex{}
	mutex.Lock()
	node.ReadCount++
	mutex.Unlock()
}

/*
   GrootGraph is the variation graph implementation used by GROOT
*/
type GrootGraph struct {
	GraphID     int
	SortedNodes []*GrootGraphNode
	NodeLookup  map[int]int
	Paths       map[int][]byte // lookup to relate PathIDs in each node to a path name
	Lengths     map[int]int    // lengths of sequences held in graph (lookup key corresponds to key in Paths)
}

// CreateGrootGraph is a GrootGraph constructor that takes a GFA instance and stores the info as a graph and then runs a topological sort
func CreateGrootGraph(gfaInstance *gfa.GFA, id int) (*GrootGraph, error) {
	// construct an empty graph
	newGraph := &GrootGraph{
		GraphID:    id,
		NodeLookup: make(map[int]int),
		Paths:      make(map[int][]byte),
		Lengths:    make(map[int]int),
	}
	// collect all the segments from the GFA instance and create the nodes
	segments, err := gfaInstance.GetSegments()
	if err != nil {
		return nil, err
	}
	for nodeIterator, segment := range segments {
		// check the segment name can be stored as an int
		id, err := strconv.Atoi(string(segment.Name))
		if err != nil {
			return nil, fmt.Errorf("could not convert segment name from GFA into an int for groot graph: %v", segment.Name)
		}
		// convert all bases to upperCase and check for non-ACTGN chars
		seq := seqio.Sequence{Seq: segment.Sequence}
		if err := seq.BaseCheck(); err != nil {
			return nil, err
		}
		newNode := &GrootGraphNode{
			SegmentID: id,
			Sequence:  seq.Seq,
			Position:  make(map[int]int),
		}
		// store the new node in the graph and record it's location in the silce by using the NodeLookup map
		newGraph.SortedNodes = append(newGraph.SortedNodes, newNode)
		newGraph.NodeLookup[id] = nodeIterator
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
			return nil, fmt.Errorf("could not convert segment name from GFA into an int for groot graph: %v", link.From)
		}
		toSegID, err := strconv.Atoi(string(link.To))
		if err != nil {
			return nil, fmt.Errorf("could not convert segment name from GFA into an int for groot graph: %v", link.To)
		}
		// add the outEdges
		nodeLocator := newGraph.NodeLookup[fromSegID]
		newGraph.SortedNodes[nodeLocator].OutEdges = append(newGraph.SortedNodes[nodeLocator].OutEdges, toSegID)
	}
	// collect all the paths from the GFA instance and add pathIDs to each node
	paths, err := gfaInstance.GetPaths()
	if err != nil {
		return nil, err
	}
	for pathIterator, path := range paths {
		// add the path name to the lookup
		newGraph.Paths[pathIterator] = path.PathName
		for _, seg := range path.SegNames {
			// strip the plus
			seg = bytes.TrimSuffix(seg, []byte("+"))
			id, err := strconv.Atoi(string(seg))
			if err != nil {
				return nil, fmt.Errorf("could not convert segment name from GFA into an int for groot graph: %v", string(seg))
			}
			nodeLocator := newGraph.NodeLookup[id]
			newGraph.SortedNodes[nodeLocator].PathIDs = append(newGraph.SortedNodes[nodeLocator].PathIDs, pathIterator)
		}
	}
	// return without toposort if only one node present (graph with single sequence)
	if len(newGraph.SortedNodes) > 1 {
		err = newGraph.topoSort()
	}
	// get and store the lengths of each sequence held in the graph
	// this function also adds Position to each node of the graph
	for pathID, _ := range newGraph.Paths {
		newGraph.Lengths[pathID] = len(newGraph.Graph2Seq(pathID))
	}
	// return the new GrootGraph
	return newGraph, err
}

// topoSort runs a topological sort on the GrootGraph
func (graph *GrootGraph) topoSort() error {
	// copy all of the graph nodes into a map (so we can keep track of what we have processed)
	nodeMap := make(map[int]*GrootGraphNode)
	toposortStart := []int{}
	seenPaths := make(map[int]int)
	for _, node := range graph.SortedNodes {
		if len(seenPaths) == len(graph.Paths) {
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
	graph.SortedNodes = []*GrootGraphNode{}
	graph.NodeLookup = make(map[int]int)
	// run the topological sort  - try starting from each node that was in the first slot of the nodeholder (start of the MSA)
	seen := make(map[int]struct{})
	for len(nodeMap) > 1 {
		for _, start := range toposortStart {
			if _, ok := nodeMap[start]; !ok {
				continue
			}
			graph.traverse(nodeMap[start], nodeMap, seen)
		}
	}
	// check all traversals have been taken
	if len(nodeMap) > 0 {
		return fmt.Errorf("topological sort failed - too many nodes remaining in the pre-sort list")
	}
	return nil
}

//  traverse is a helper method to topologically sort a graph
func (graph *GrootGraph) traverse(node *GrootGraphNode, nodeMap map[int]*GrootGraphNode, seen map[int]struct{}) {
	// skip if we are already handling the current node
	if _, ok := seen[node.SegmentID]; ok {
		return
	}
	// make sure node is still in the graph
	if _, ok := nodeMap[node.SegmentID]; ok {
		// record that we are processing this node
		seen[node.SegmentID] = struct{}{}
		// sort the output nodes for this node and then traverse them in reverse order
		sort.Sort(sort.Reverse(sort.IntSlice(node.OutEdges)))
		for i, j := 0, len(node.OutEdges); i < j; i++ {
			// check if the outedges have been traversed
			if _, ok := nodeMap[node.OutEdges[i]]; !ok {
				continue
			}
			graph.traverse(nodeMap[node.OutEdges[i]], nodeMap, seen)
		}
		// delete the node from the temporary holders
		delete(nodeMap, node.SegmentID)
		delete(seen, node.SegmentID)
		// update the sorted node slice and the lookup
		graph.SortedNodes = append([]*GrootGraphNode{node}, graph.SortedNodes...)
		graph.NodeLookup[node.SegmentID] = len(nodeMap)
	}
}

// Graph2Seq is a method to convert a variation graph to linear reference sequences
func (graph *GrootGraph) Graph2Seq(pathID int) []byte {
	newSequence := []byte{}
	for _, node := range graph.SortedNodes {
		for _, id := range node.PathIDs {
			if id == pathID {
				// update the Position field of the current node
				node.Position[pathID] = len(newSequence)
				// build the sequence
				newSequence = append(newSequence, node.Sequence...)
			}
		}
	}
	return newSequence
}

// WindowGraph is a method to window every node of the graph, generate minhash signautes and return unique ones and their graph locations
func (graph *GrootGraph) WindowGraph(windowSize, kSize, sigSize int) <-chan *Window {
	// generate a hashing function to convert signatures (uint64[]) to a string
	sigHasher := GenSigHasher(2)
	// make the channel to send windows over
	windowChan := make(chan *Window)
	go func() {
		defer close(windowChan)
		// for each reference in the graph, get all the windows possible from every base in every node
		for pathID := range graph.Paths {
			for _, node := range graph.SortedNodes {
				for offset := 0; offset < len(node.Sequence); offset++ {
					var wg sync.WaitGroup
					sendPath := make(chan []byte)
					// launch the DFS
					wg.Add(1)
					go func() {
						defer wg.Done()
						_ = graph.recursiveDFS(node, offset, []byte{}, sendPath, windowSize, pathID)
					}()
					go func() {
						wg.Wait()
						close(sendPath)
					}()
					// sigChecker sees if we get multiple signatures for the same node+offset
					sigChecker := make(map[string]struct{})
					// get MinHash signature for all the windows possible from the current base in the current node
					for path := range sendPath {
						// get a minhash signature
						windowSeq := seqio.Sequence{Seq: path}
						sig := windowSeq.RunMinHash(kSize, sigSize).Signature()
						// check if we have seen this signature for this node and offset
						hashedSig := sigHasher(sig)
						if _, ok := sigChecker[hashedSig]; !ok {
							// add it to the checker
							sigChecker[hashedSig] = struct{}{}
							// create a GrootGraph window
							newWindow := &Window{
								GraphID: graph.GraphID,
								Node:    node.SegmentID,
								OffSet:  offset,
								Sig:     sig,
							}
							// send window
							windowChan <- newWindow
						}
					}
				}
			}
		}
	}()
	return windowChan
}

func (graph *GrootGraph) recursiveDFS(currentNode *GrootGraphNode, offset int, path []byte, sendPath chan []byte, windowSize, pathID int) bool {
	// a switch to halt the current path traversal upon exact alignment covering whole read
	windowComplete := false
	// end the DFS if this node is not derived from the pathID
	nodeOK := false
	for _, ID := range currentNode.PathIDs {
		if ID == pathID {
			nodeOK = true
			break
		}
	}
	if nodeOK == false {
		return false
	}
	// iterate over the sequence
	for _, base := range currentNode.Sequence[offset:] {
		// update the path
		path = append(path, base)
		// end the DFS if we have traversed enough bases to equal the window size
		if len(path) == windowSize {
			pathCopy := make([]byte, len(path))
			for i, j := range path {
				pathCopy[i] = j
			}
			sendPath <- pathCopy
			return true
		}
	}
	// end the DFS if there are no out edges
	if len(currentNode.OutEdges) == 0 {
		return false
	}
	// else, continue along the current traversal
	for _, neighbourNode := range currentNode.OutEdges {
		nodeLocator := graph.NodeLookup[neighbourNode]
		nextNode := graph.SortedNodes[nodeLocator]
		if result := graph.recursiveDFS(nextNode, 0, path, sendPath, windowSize, pathID); result == true {
			windowComplete = true
		}
	}
	return windowComplete
}

/// GenSigHasher generates a hashing function to convert a MinHash signature to a string so that we can use it as a map key and check for duplicates quickly
func GenSigHasher(hashValueSize int) func([]uint64) string {
	return func(sig []uint64) string {
		hashedSig := make([]byte, hashValueSize*len(sig))
		buf := make([]byte, 8)
		for i, v := range sig {
			// use the ByteOrder interface to write binary data
			// use the LittleEndian implementation and call the Put method
			binary.LittleEndian.PutUint64(buf, v)
			copy(hashedSig[i*hashValueSize:(i+1)*hashValueSize], buf[:hashValueSize])
		}
		return string(hashedSig)
	}
}

// DumpGraph is a method to save a GrootGraph in GFA format
func (graph *GrootGraph) DumpGraph(dirName string) (int, error) {
	// a flag to prevent dumping graphs which had no reads align
	graphUsed := false
	t := time.Now()
	stamp := fmt.Sprintf("variation graph created by groot at: %v", t.Format("Mon Jan _2 15:04:05 2006"))
	// create a GFA instance
	newGFA := gfa.NewGFA()
	_ = newGFA.AddVersion(1)
	newGFA.AddComment([]byte(stamp))
	// transfer all the GrootGraphNode content to the GFA instance
	for _, node := range graph.SortedNodes {
		// record if this graph has had reads align
		if (graphUsed == false) && (node.ReadCount > 0) {
			graphUsed = true
		}
		segID := []byte(strconv.Itoa(node.SegmentID))
		// create the segment
		seg, err := gfa.NewSegment(segID, []byte(node.Sequence))
		if err != nil {
			return 0, err
		}
		// I'm weighting the RC so that bandage will display the depth nicely for the graph paths
		readCount := fmt.Sprintf("RC:i:%d", (node.ReadCount * len(node.Sequence)))
		// FragCount is the raw readcount
		fragCount := fmt.Sprintf("FC:i:%d", node.ReadCount)
		ofs, err := gfa.NewOptionalFields([]byte(readCount), []byte(fragCount))
		if err != nil {
			return 0, err
		}
		seg.AddOptionalFields(ofs)
		seg.Add(newGFA)
		// create the links
		for _, outEdge := range node.OutEdges {
			toSeg := []byte(strconv.Itoa(outEdge))
			link, err := gfa.NewLink(segID, []byte("+"), toSeg, []byte("+"), []byte("0M"))
			if err != nil {
				return 0, err
			}
			link.Add(newGFA)
		}
	}
	// don't save the graph if no reads aligned
	if graphUsed == false {
		return 0, nil
	}
	// create the paths
	for pathID, pathName := range graph.Paths {
		segments, overlaps := [][]byte{}, [][]byte{}
		for _, node := range graph.SortedNodes {
			for _, id := range node.PathIDs {
				if id == pathID {
					segment := strconv.Itoa(node.SegmentID) + "+"
					overlap := strconv.Itoa(len(node.Sequence)) + "M"
					segments = append(segments, []byte(segment))
					overlaps = append(overlaps, []byte(overlap))
					break
				}
			}
		}
		// add the path
		path, err := gfa.NewPath(pathName, segments, overlaps)
		if err != nil {
			return 0, err
		}
		path.Add(newGFA)
	}
	// create a gfaWriter and write the GFA instance
	firstSeqID := graph.Paths[0]
	if bytes.Contains(firstSeqID, []byte("/")) {
		firstSeqID = bytes.Replace(firstSeqID, []byte("/"), []byte(""), -1)
	}
	if bytes.Contains(firstSeqID, []byte(":")) {
		firstSeqID = bytes.Replace(firstSeqID, []byte(":"), []byte(""), -1)
	}
	fileName := fmt.Sprintf("%v/%v-groot-graph.gfa", dirName, string(firstSeqID))
	outfile, err := os.Create(fileName)
	if err != nil {
		return 0, err
	}
	defer outfile.Close()
	writer, err := gfa.NewWriter(outfile, newGFA)
	if err != nil {
		return 0, err
	}
	err = newGFA.WriteGFAContent(writer)
	return 1, nil
}

/*
  A struct to store multiple graphs
*/
type GraphStore map[int]*GrootGraph

// Dump is a method to save a GrootGraph to file
func (graphStore *GraphStore) Dump(path string) error {
	file, err := os.Create(path)
	if err == nil {
		encoder := gob.NewEncoder(file)
		encoder.Encode(graphStore)
	}
	file.Close()
	return err
}

// Load is a method to load a GrootGraph from file
func (graphStore *GraphStore) Load(path string) error {
	file, err := os.Open(path)
	if err == nil {
		decoder := gob.NewDecoder(file)
		err = decoder.Decode(graphStore)
	}
	file.Close()
	return err
}

// GetRefs is a method to convert all paths held in graphStore to sam.References
func (graphStore GraphStore) GetRefs() (map[int][]*sam.Reference, error) {
	references := make(map[int][]*sam.Reference)
	for graphID, grootGraph := range graphStore {
		references[graphID] = make([]*sam.Reference, len(grootGraph.Paths))
		for pathID, path := range grootGraph.Paths {
			reference, err := sam.NewReference(string(path), "", "", grootGraph.Lengths[pathID], nil, nil)
			if err != nil {
				return nil, err
			}
			references[graphID][pathID] = reference
		}
	}
	return references, nil
}

/*
  A struct to hold graph window information (for minhashing)
*/
type Window struct {
	GraphID int
	Node    int
	OffSet  int
	Sig     []uint64
}
