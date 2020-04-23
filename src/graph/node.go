package graph

// Nodes is a slice of GrootGraphNodes
type Nodes []uint64

// methods for Nodes to satisfy the sort interface
func (Nodes Nodes) Len() int           { return len(Nodes) }
func (Nodes Nodes) Swap(i, j int)      { Nodes[i], Nodes[j] = Nodes[j], Nodes[i] }
func (Nodes Nodes) Less(i, j int) bool { return Nodes[i] < Nodes[j] }

// GrootGraphNode is a GFA segment (plus the extra info from path, links etc.)
// Note: GrootGraphNodes are not set up for concurrent access
type GrootGraphNode struct {
	SegmentID     uint64
	SegmentLength float64
	Sequence      []byte
	OutEdges      Nodes
	PathIDs       []uint32    // PathIDs are the lookup IDs to the linear reference sequences that use this segment (value corresponds to key in GrootGraph.Paths)
	Position      map[int]int // the start position of this segment in each reference sequence, using 1-based indexing (lookup key corresponds to PathID)
	KmerFreq      float64
	Marked        bool // TODO: tmp idea to mark nodes during pruning, rather than deleting
}

// IncrementKmerFreq is a method to increment a node's k-mer count
func (GrootGraphNode *GrootGraphNode) IncrementKmerFreq(increment float64) error {
	GrootGraphNode.KmerFreq += increment
	return nil
}
