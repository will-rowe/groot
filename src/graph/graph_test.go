/*
	tests for the graph package
*/
package graph

import (
	"fmt"
	"testing"
)

// test intput
var (
	seq1 = []byte("---ACT")
	seq2 = []byte("---AGT")
	seq3 = []byte("CCCACT")
	seq4 = []byte("G-CACT")
	seq5 = []byte("-ACT-GC--")
	seq6 = []byte("AACTAGC--")
)

// expected edges for the sort graph from test function 2
var testInEdges = [][]int{
	{},
	{0},
	{0},
	{1, 4},
}
var testOutEdges = [][]int{
	{4, 1},
	{2},
	{2},
	{},
}

// the order of nodes expected in the sorted graph from test function 3
var testNodeOrder = [11]int{12, 6, 7, 8, 0, 1, 4, 2}

// test function to check equality of slices
func byteSliceCheck(a, b []byte) bool {
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// this test makes sure that the graph can be initialised using the representative sequence (seq1)
func TestBasicGraph(t *testing.T) {
	// create the graph
	graph, err := NewGraph(seq1)
	if err != nil {
		t.Fatalf("could not create the graph: %v\n", err)
	}
	// topologically sort it
	if err := graph.TopSort(); err != nil {
		t.Fatalf("could not run topological sort: %v\n", err)
	}
	// run Graph2Seq and make sure it matches the input sequence (ignore leading -)
	testSeq, _, err := graph.Graph2Seq(0)
	if err != nil {
		t.Fatalf("could not extract sequence from the graph: %v\n", err)
	}
	if byteSliceCheck(seq1[3:], testSeq) != true {
		t.Fatalf("graph2Seq method did not reproduce input representative sequence\na: %v\nb: %v\n", string(seq1[3:]), string(testSeq))
	}
}

// this test makes sure that the graph can have variant nodes added using the MSA
func TestVariantGraph(t *testing.T) {
	// create the graph
	graph, err := NewGraph(seq1)
	if err != nil {
		t.Fatalf("could not create the graph: %v\n", err)
	}
	// add in a second sequence from the MSA
	if err = graph.AddVariantNodes(seq2, 1); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	fmt.Println("variant graph nodes:")
	fmt.Println(graph.NodeHolder)
	// check the topological sort
	if err := graph.TopSort(); err != nil {
		t.Fatalf("could not run topological sort: %v\n", err)
	}
	fmt.Println("sorted variant graph:")
	fmt.Println(graph.SortedNodes)
	// check all the edges are correct
	for i, nodeEdges := range testInEdges {
		for j, edge := range nodeEdges {
			if graph.SortedNodes[i].InEdges[j] != edge {
				t.Fatalf("incorrect IN edges\n")
			}
		}
	}
	for i, nodeEdges := range testOutEdges {
		for j, edge := range nodeEdges {
			if graph.SortedNodes[i].OutEdges[j] != edge {
				t.Fatalf("incorrect OUT edges\n")
			}
		}
	}
	// check for nodes correctly labelled as shared
	if len(graph.SortedNodes[0].Parent) != 2 && len(graph.SortedNodes[3].Parent) != 2 {
		t.Fatalf("The first and last nodes of the sorted should be marked as shared\n")
	}
}

// this test makes sure that the graph can have variant nodes added using the MSA
func TestMoreComplexVariantGraph(t *testing.T) {
	// create the graph
	graph, err := NewGraph(seq1)
	if err != nil {
		t.Fatalf("could not create the graph: %v\n", err)
	}
	// add in additional sequences from the MSA
	if err = graph.AddVariantNodes(seq2, 1); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	if err = graph.AddVariantNodes(seq3, 2); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	if err = graph.AddVariantNodes(seq4, 3); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	fmt.Println("variant graph nodes:")
	for i, slot := range graph.NodeHolder {
		fmt.Println(i, slot)
	}
	// check the topological sort
	if err := graph.TopSort(); err != nil {
		t.Fatalf("could not run topological sort: %v\n", err)
	}
	fmt.Println("sorted variant graph:")
	for i, node := range graph.SortedNodes {
		fmt.Printf("Node ID: %d\tParent(s): %d\tBase: %v\tIn edges: %d\tOut edges: %d\n", node.ID, node.Parent, string(node.Base), node.InEdges, node.OutEdges)
		if node.ID != testNodeOrder[i] {
			t.Fatalf("topological sort is incorrect")
		}
	}
	// check the representative sequence can be generated
	// run Graph2Seq and make sure it matches the input sequence (ignore leading -)
	testSeq, _, err := graph.Graph2Seq(0)
	if err != nil {
		t.Fatalf("could not return sequence from the graph: %v\n", err)
	}
	if byteSliceCheck(seq1[3:], testSeq) != true {
		t.Fatalf("graph2Seq method did not reproduce input representative sequence\na: %v\nb: %v\n", string(seq1[3:]), string(testSeq))
	}
}

func TestGraph2Seq(t *testing.T) {
	// create the graph
	graph, err := NewGraph(seq1)
	if err != nil {
		t.Fatalf("could not create the graph: %v\n", err)
	}
	// add in additional sequences from the MSA
	if err = graph.AddVariantNodes(seq2, 1); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	if err = graph.AddVariantNodes(seq3, 2); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	if err = graph.AddVariantNodes(seq4, 3); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	// check the topological sort
	if err := graph.TopSort(); err != nil {
		t.Fatalf("could not run topological sort: %v\n", err)
	}
	// check each reference sequence can be extracted from the graph
	testSeq1, startNode, err := graph.Graph2Seq(0)
	if err != nil {
		t.Fatalf("could not extract seq1 from the graph: %v\n", err)
	}
	if byteSliceCheck(seq1[3:], testSeq1) != true {
		t.Fatalf("graph2Seq method did not reproduce seq1\na: %v\nb: %v\n", string(seq1[3:]), string(testSeq1))
	}
	if startNode != 0 {
		t.Fatalf("sequence extracted but incorrect start node returned")
	}
	testSeq2, startNode, err := graph.Graph2Seq(1)
	if err != nil {
		t.Fatalf("could not extract seq2 from the graph: %v\n", err)
	}
	if byteSliceCheck(seq2[3:], testSeq2) != true {
		t.Fatalf("graph2Seq method did not reproduce seq2\na: %v\nb: %v\n", string(seq2[3:]), string(testSeq2))
	}
	if startNode != 0 {
		t.Fatalf("sequence extracted but incorrect start node returned")
	}
	testSeq3, startNode, err := graph.Graph2Seq(2)
	if err != nil {
		t.Fatalf("could not extract seq3 from the graph: %v\n", err)
	}
	if byteSliceCheck(seq3, testSeq3) != true {
		t.Fatalf("graph2Seq method did not reproduce seq3\na: %v\nb: %v\n", string(seq3), string(testSeq3))
	}
	if startNode != 6 {
		t.Fatalf("sequence extracted but incorrect start node returned")
	}
	// check the window parameter works
	testSeq3window, _, err := graph.Graph2Seq(2, 3, 3)
	if err != nil {
		t.Fatalf("could not extract window of seq3 from the graph: %v\n", err)
	}
	if byteSliceCheck(seq3[3:], testSeq3window) != true {
		t.Fatalf("graph2Seq method did not reproduce window of seq3\na: %v\nb: %v\n", string(seq3[3:]), string(testSeq3window))
	}
}

// this test tries a gapped representative sequence
func TestGappedRepresentativeSeq(t *testing.T) {
	// create the graph
	graph, err := NewGraph(seq5)
	if err != nil {
		t.Fatalf("could not create the graph: %v\n", err)
	}
	// add in additional sequences from the MSA
	if err = graph.AddVariantNodes(seq6, 1); err != nil {
		t.Fatalf("could not add variants the graph: %v\n", err)
	}
	// check the topological sort
	if err := graph.TopSort(); err != nil {
		t.Fatalf("could not run topological sort: %v\n", err)
	}
	fmt.Println("sorted variant graph:")
	for _, node := range graph.SortedNodes {
		fmt.Printf("Node ID: %d\tParent(s): %d\tBase: %v\tIn edges: %d\tOut edges: %d\n", node.ID, node.Parent, string(node.Base), node.InEdges, node.OutEdges)
	}
}
