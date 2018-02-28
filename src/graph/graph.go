// the graph package stores a reference sequence and its variants as a directed multigraph
package graph

import (
	"encoding/gob"
	"errors"
	"os"
	"sort"

	"github.com/biogo/hts/sam"
)

/*
  Structs to store reference sequence (+variants) as a graph
*/
type Node struct {
	ID       int
	Parent   []int
	Base     byte
	InEdges  []int
	OutEdges []int
	Position map[int]int
}
type Graph struct {
	NodeHolder  [][]Node
	NodeTotal   int
	SortedNodes []Node
	NodeLookUp  map[int]int
	ARGs        [][]byte
	Lengths     []int
}

// method to count the number of nodes in a map and store the value
func (self *Graph) countNodes() {
	self.NodeTotal = 0
	for _, nodes := range self.NodeHolder {
		for range nodes {
			self.NodeTotal++
		}
	}
}

// method to add additional sequences to a graph
func (self *Graph) AddVariantNodes(refSeq []byte, refID int) error {
	if len(refSeq) != len(self.NodeHolder) {
		return errors.New("can't add variant nodes - sequence is a different length to the sequence used to initialise the graph")
	}
	// refresh the current graph node count - use this to create ID for new nodes
	self.countNodes()
	// seqPos is the position of the node in the sequence (not the MSA) - it uses 0-base index and is used to report the alignment location relative to the sequence, not the graph
	seqPos := 0
	finalNodeSlot := 0
	// iterate over the positions in the MSA -
	for msaPos, base := range refSeq {
		// skip deletions
		if base == 45 {
			continue
		}
		newNodeID := self.NodeTotal
		// connect the node to the nodes from the rest of this sequence
		varIn, varOut := []int{}, []int{}
		if seqPos != 0 {
			varIn = append(varIn, newNodeID-1)
		}
		varOut = append(varOut, newNodeID+1)
		// create the new node
		varNode := NewNode(newNodeID, refID, base, varIn, varOut)
		// add the sequence position to the node
		varNode.Position = make(map[int]int)
		varNode.Position[refID] = seqPos
		seqPos++
		// add the node to the graph node holder - use the slot corresponding to the current MSA position
		self.NodeHolder[msaPos] = append(self.NodeHolder[msaPos], varNode)
		// increment the node tally for the graph
		self.NodeTotal++
		finalNodeSlot = msaPos
	}
	// graph clean up - remove OutEdges of the final node
	self.NodeHolder[finalNodeSlot][len(self.NodeHolder[finalNodeSlot])-1].OutEdges = nil
	return nil
}

// TODO: this needs some checking adding
// method to prune duplicated nodes and re-order edges
func (self *Graph) Prune() error {
	// tmp store of consolidated nodes
	rmNodes := make(map[int]int)
	// iterate over the node holder
	for slot := 0; slot < len(self.NodeHolder); slot++ {
		// skip empty slots
		if len(self.NodeHolder[slot]) == 0 {
			continue
		}
		// iterate over the nodes in this slot, identifying any duplicate nodes (same base)
		duplicates := make(map[string]int)
		for nodeLocator, currentNode := range self.NodeHolder[slot] {
			if _, ok := duplicates[string(currentNode.Base)]; !ok {
				duplicates[string(currentNode.Base)] = nodeLocator
			} else {
				// the current node is a duplicate - the previous occurrence is now the node we will consolidate the duplicate(s) to
				consolidator := duplicates[string(currentNode.Base)]
				// update the consolidator to contain the current nodes in/out edges + parent identifier
				self.NodeHolder[slot][consolidator].InEdges = append(self.NodeHolder[slot][consolidator].InEdges, currentNode.InEdges...)
				self.NodeHolder[slot][consolidator].OutEdges = append(self.NodeHolder[slot][consolidator].OutEdges, currentNode.OutEdges...)
				self.NodeHolder[slot][consolidator].Parent = append(self.NodeHolder[slot][consolidator].Parent, currentNode.Parent...)

				// add all the start position of the current node to the consolidation node's positions map (providing they aren't gap nodes)
				self.NodeHolder[slot][consolidator].Position[currentNode.Parent[0]] = currentNode.Position[currentNode.Parent[0]]
				// finally, remove the duplicated node from the slot in the holder (blank it for now, toposort will then remove later) and update the rmNodes map
				self.NodeHolder[slot][nodeLocator] = Node{ID: -1}
				rmNodes[currentNode.ID] = self.NodeHolder[slot][consolidator].ID
			}
		}
	}
	// iterate over the node holder again, tidy up all references to removed nodes and remove duplicate edges
	for slot := 0; slot < len(self.NodeHolder); slot++ {
		for nodeLocator := range self.NodeHolder[slot] {
			// update in edges + remove duplicates
			inDups := make(map[int]struct{})
			updatedInEdges := []int{}
			for i := 0; i < len(self.NodeHolder[slot][nodeLocator].InEdges); i++ {
				inEdge := self.NodeHolder[slot][nodeLocator].InEdges[i]
				if consolidator, ok := rmNodes[inEdge]; ok {
					self.NodeHolder[slot][nodeLocator].InEdges[i] = consolidator
					inEdge = consolidator
				}
				// if this is a duplicate node ID,  flag it
				if _, ok := inDups[inEdge]; !ok {
					inDups[inEdge] = struct{}{}
					updatedInEdges = append(updatedInEdges, inEdge)
				}
			}
			self.NodeHolder[slot][nodeLocator].InEdges = updatedInEdges
			// update out edges + remove duplicates
			outDups := make(map[int]struct{})
			updatedOutEdges := []int{}
			for i := 0; i < len(self.NodeHolder[slot][nodeLocator].OutEdges); i++ {
				outEdge := self.NodeHolder[slot][nodeLocator].OutEdges[i]
				if consolidator, ok := rmNodes[outEdge]; ok {
					self.NodeHolder[slot][nodeLocator].OutEdges[i] = consolidator
					outEdge = consolidator
				}
				// if this is a duplicate node ID, flag it
				if _, ok := outDups[outEdge]; !ok {
					outDups[outEdge] = struct{}{}
					updatedOutEdges = append(updatedOutEdges, outEdge)
				}
			}
			self.NodeHolder[slot][nodeLocator].OutEdges = updatedOutEdges
		}
	}
	return nil
}

// method to topologically sort the graph
func (self *Graph) TopSort() error {
	// currently topsort can only be run once, after adding all variant nodes
	if len(self.SortedNodes) != 0 {
		return errors.New("this graph has already been topologically sorted")
	}
	if err := self.Prune(); err != nil {
		return errors.New("could not prune the graph")
	}
	self.countNodes()
	// copy all of the graph nodes into a map (so we can keep track of what we have processed)
	nodeMap := make(map[int]Node)
	toposortStart := []int{}
	seenParent := make(map[int]int)
	for _, nodes := range self.NodeHolder {
		// skip the current nodeholder position if it contains no nodes
		if len(nodes) == 0 {
			continue
		}
		for _, node := range nodes {
			// ignore empty nodes (created during pruning)
			if node.ID == -1 {
				self.NodeTotal--
				continue
			}
			// if this node is from a sequence we have not seen yet, keep the node.ID as a starting node for toposort
			for _, parent := range node.Parent {
				if _, ok := seenParent[parent]; !ok {
					toposortStart = append(toposortStart, node.ID)
				}
			}
			// check for duplicate nodes
			if _, ok := nodeMap[node.ID]; ok {
				return errors.New("graph contains duplicate nodes")
			}
			nodeMap[node.ID] = node
		}
	}

	// make a lookup to relate the ordering of the sorted nodes to the node ID
	self.NodeLookUp = make(map[int]int, self.NodeTotal)
	// check we have added all the nodes to the map and then delete the node holder
	if len(nodeMap) != self.NodeTotal {
		return errors.New("node count error during TopSort")
	} else {
		self.NodeHolder = nil
	}
	// run the topological sort  - try starting from each node that was in the first slot of the nodeholder (start of the MSA)
	seen := make(map[int]struct{})
	for len(nodeMap) > 1 {
		for _, start := range toposortStart {
			if _, ok := nodeMap[start]; !ok {
				continue
			}
			traverse(nodeMap[start], nodeMap, seen, self)
		}
	}
	// check all traversals have been taken
	if len(nodeMap) > 1 {
		return errors.New("topological sort failed - too many nodes remaining in the pre-sort list")
	}
	// add the final node to the sorted slice, sorting the out edges
	for _, node := range nodeMap {
		sort.Sort(sort.Reverse(sort.IntSlice(node.OutEdges)))
		self.SortedNodes = append([]Node{node}, self.SortedNodes...)
		self.NodeLookUp[node.ID] = 0
	}
	return nil
}

//  helper function to topologically sort a graph
func traverse(node Node, nodeMap map[int]Node, seen map[int]struct{}, graph *Graph) {
	// skip if we are already handling the current node
	if _, ok := seen[node.ID]; ok {
		return
	}
	// make sure node is still in the graph
	if _, ok := nodeMap[node.ID]; ok {
		// record that we are processing this node
		seen[node.ID] = struct{}{}
		// sort the output nodes for this node and then traverse them in reverse order
		sort.Sort(sort.Reverse(sort.IntSlice(node.OutEdges)))
		for i, j := 0, len(node.OutEdges); i < j; i++ {
			// check if the outedges have been traversed
			if _, ok := nodeMap[node.OutEdges[i]]; !ok {
				continue
			}
			traverse(nodeMap[node.OutEdges[i]], nodeMap, seen, graph)
		}

		// delete the node from the temporary holders
		delete(nodeMap, node.ID)
		delete(seen, node.ID)
		// update the sorted node slice and the lookup
		graph.SortedNodes = append([]Node{node}, graph.SortedNodes...)
		graph.NodeLookUp[node.ID] = len(nodeMap)
	}
}

// method to send the a reference held in graph to a sequence (as []byte)
// it also returns the node ID of the first base in the returned sequence
// requires the referenceID of the sequence to extract from the graph
// optional: two integers to specify a window to extract (the start position in the sequence (0 based) and the length of the window)
func (self *Graph) Graph2Seq(referenceID int, window ...int) ([]byte, int, error) {
	var seq []byte
	windowMode := false
	nodeCounter, startNode := 0, 0
	// make sure the graph has been sorted
	if len(self.SortedNodes) == 0 {
		return nil, 0, errors.New("the graph needs to be topologically sorted first")
	}
	if len(window) > 2 {
		return nil, 0, errors.New("Graph2Seq only accepts a maximum of 2 integers for the window parameter (start position and window length)")
	}
	if len(window) == 2 {
		windowMode = true
	}
	// loop through the sorted nodes of the graph
	for i, j := 0, len(self.SortedNodes); i < j; i++ {
		// for each node, loop through the parent sequences
		for k, l := 0, len(self.SortedNodes[i].Parent); k < l; k++ {
			// check the current node if the parent sequence matches the referenceID of the sequence we are trying to build
			if self.SortedNodes[i].Parent[k] == referenceID {
				// increment the nodeCounter
				nodeCounter++
				// skip the node if we are windowing and this node is not in the sequence window
				if windowMode == true && nodeCounter-1 < window[0] {
					continue
				}
				// add the base to the sequence
				seq = append(seq, self.SortedNodes[i].Base)
				// record the starting node of the sequence
				if len(seq) == 1 {
					startNode = self.SortedNodes[i].ID
				}
				// if windowing and the seq is now the right length, return it
				if windowMode == true && len(seq) == window[1] {
					return seq, startNode, nil
				}
			}
		}
	}
	// check we have extracted something
	if len(seq) == 0 {
		return nil, 0, errors.New("could not extract reference sequence from graph - check the correct ID was provided")
	}
	// if windowing, check aren't about to send back something the wrong length (right length should have already been sent back)
	if windowMode == true && len(seq) != window[1] {
		return nil, 0, errors.New("could not extract a window of the specified length - ran out of graph")
	}
	return seq, startNode, nil
}

/*
  A function to construct a new node
*/
func NewNode(id int, parent int, base byte, in []int, out []int) Node {
	node := new(Node)
	node.ID = id
	node.Parent = []int{parent}
	node.Base = base
	node.InEdges = in
	node.OutEdges = out
	return *node
}

/*
  A function to construct an empty graph
*/
func InitGraph(nhLength int) Graph {
	graph := new(Graph)
	nodeHolder := make([][]Node, nhLength)
	graph.NodeHolder = nodeHolder
	graph.SortedNodes = []Node{}
	return *graph
}

/*
  A function to add nodes to a graph from a representative sequence
*/
func NewGraph(refSeq []byte) (*Graph, error) {
	// create an empty graph
	graph := InitGraph(len(refSeq))
	// nodeID is used as a counter that equates to how many nodes are currently in the graph
	nodeID := 0
	// seqPos is the position of the node in the sequence (not the MSA) - it uses 0-base index and is used to report the alignment location relative to the sequence, not the graph
	seqPos := 0
	// finalNode records which slot in the NodeHolder the final node was added for this seq
	finalNode := 0
	// create the nodes by iterating over the length of the MSA
	for msaPos, base := range refSeq {
		// skip deletions
		if base == 45 {
			continue
		}
		// if current node is not the first in the graph, add the in edges
		in, out := []int{}, []int{}
		if nodeID != 0 {
			in = append(in, nodeID-1)
		}
		// add the out edges (if it is the last node added, we need to remove the out edge at the end of the new graph function)
		out = append(out, nodeID+1)
		// create a node
		node := NewNode(nodeID, 0, base, in, out)
		// add the sequence position to the node
		node.Position = make(map[int]int)
		node.Position[0] = seqPos
		seqPos++
		// add the node to the graph node holder (graph prior to topological sort)
		graph.NodeHolder[msaPos] = append(graph.NodeHolder[msaPos], node)
		// increment the nodeID
		nodeID++
		finalNode = msaPos
	}
	// remove the out edge from the final node added during the initial graph build
	graph.NodeHolder[finalNode][0].OutEdges = nil
	// count the number of nodes, check the graph and return it
	graph.countNodes()
	if graph.NodeTotal == 0 {
		return nil, errors.New("could create any nodes from the input sequence")
	}
	return &graph, nil
}

/*
  A struct to store multiple graphs
*/
type GraphStore map[int]*Graph

// method to dump the graph to file
func (self *GraphStore) Dump(path string) error {
	file, err := os.Create(path)
	if err == nil {
		encoder := gob.NewEncoder(file)
		encoder.Encode(self)
	}
	file.Close()
	return err
}

// method to load a graph from file into a graph struct
func (self *GraphStore) Load(path string) error {
	file, err := os.Open(path)
	if err == nil {
		decoder := gob.NewDecoder(file)
		err = decoder.Decode(self)
	}
	file.Close()
	return err
}

// method to convert all ARGs held in graphStore to sam.References
func (self GraphStore) GetRefs() (map[int][]*sam.Reference, error) {
	references := make(map[int][]*sam.Reference, len(self))
	for graphIterator, j := 0, len(self); graphIterator < j; graphIterator++ {
		references[graphIterator] = make([]*sam.Reference, len(self[graphIterator].ARGs))
		for refIterator, l := 0, len(self[graphIterator].ARGs); refIterator < l; refIterator++ {
			reference, err := sam.NewReference(string(self[graphIterator].ARGs[refIterator]), "", "", self[graphIterator].Lengths[refIterator], nil, nil)
			if err != nil {
				return nil, err
			}
			references[graphIterator][refIterator] = reference
		}
		self[graphIterator].ARGs, self[graphIterator].Lengths = nil, nil
	}
	return references, nil
}

/*
  A struct to hold graph window information (for minhashing)
*/
type Window struct {
	Graph int
	Node  int
	Sig   []uint64
}
