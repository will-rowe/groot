/*
	tests for the alignment package

    the alignment package performs the following hierarchical alignment steps:
    1. exact alignment
    2. seed shuffle + exact alignment
    3. gapped end alignment
    4. inexact alignment
*/
package alignment

import (
	"errors"
	"fmt"
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/seqio"
)

/*
  setup the tests
*/
var (
	// graph
	seq1 = []byte("AGTGGTCATA")
	seq2 = []byte("A-GGGTCA-A")
	// reads
	l1 = []byte("@read_A") // single exact match
	l2 = []byte("agtg")
	l3 = []byte("@read_B") // multiple exact match
	l4 = []byte("gtca")
	l5 = []byte("@read_D") // single inexact match
	l6 = []byte("agtggtcatt")
)

// helper funcs
func setupGraph(seq1, seq2 []byte) (graph.Graph, []*sam.Reference, error) {
	testGraph, err := graph.NewGraph(seq1)
	if err != nil {
		return graph.Graph{}, nil, errors.New(fmt.Sprintf("could not create the graph: %v\n", err))
	}
	if err = testGraph.AddVariantNodes(seq2, 1); err != nil {
		return graph.Graph{}, nil, errors.New(fmt.Sprintf("could not add variants the graph: %v\n", err))
	}
	if err := testGraph.TopSort(); err != nil {
		return graph.Graph{}, nil, errors.New(fmt.Sprintf("could not run topological sort: %v\n", err))
	}
	ref1, err := sam.NewReference("reference_A", "", "", len(seq1), nil, nil)
	if err != nil {
		return graph.Graph{}, nil, errors.New(fmt.Sprintf("could not create SAM references: %v\n", err))
	}
	ref2, err := sam.NewReference("reference_B", "", "", len(seq2), nil, nil)
	if err != nil {
		return graph.Graph{}, nil, errors.New(fmt.Sprintf("could not create SAM references: %v\n", err))
	}
	testRefs := []*sam.Reference{ref1, ref2}
	return *testGraph, testRefs, nil
}

func setupRead(la, lb, lc, ld []byte, node int) (seqio.FASTQread, error) {
	testRead, err := seqio.NewFASTQread(la, lb, lc, ld)
	if err != nil {
		return seqio.FASTQread{}, errors.New(fmt.Sprintf("could not generate FASTQ read: %v\n", err))
	}
	err = testRead.BaseCheck()
	if err != nil {
		return seqio.FASTQread{}, errors.New(fmt.Sprintf("could not generate FASTQ read: %v\n", err))
	}
	// this seed is for graph ID 0 (not needed here), node 0 (the first node in the reference graph), 100% JS and not reverse complemented
	seed := seqio.Key{0, node, false}
	testRead.Seeds = []seqio.Key{seed}
	return testRead, nil
}

/*
  begin the tests
*/
// this tests the first alignment method (exact alignment)
func TestExactMatch(t *testing.T) {
	// create the graph
	testGraph, testRefs, err := setupGraph(seq1, seq2)
	if err != nil {
		t.Fatal(err)
	}
	// create the read + seed
	testRead, err := setupRead(l1, l2, []byte("+"), []byte("===="), 0)
	if err != nil {
		t.Fatal(err)
	}
	// test the alignment
	alignments := Align(testRead, 0, &testGraph, testRefs, 5)
	collector := []*sam.Record{}
	for alignment := range alignments {
		fmt.Println(alignment)
		collector = append(collector, alignment)
	}
	if len(collector) > 1 {
		t.Fatalf("there should not be more than one alignment")
	}
	if collector[0].Ref.Name() != "reference_A" {
		t.Fatalf("aligned to wrong reference")
	}
}

// this tests the first alignment method (exact alignment) but with 2 possible matches
func TestMultipleExactMatch(t *testing.T) {
	// create the graph
	testGraph, testRefs, err := setupGraph(seq1, seq2)
	if err != nil {
		t.Fatal(err)
	}
	// create the read + seed
	testRead, err := setupRead(l3, l4, []byte("+"), []byte("===="), 4)
	if err != nil {
		t.Fatal(err)
	}
	// test the alignment
	alignments := Align(testRead, 0, &testGraph, testRefs, 5)
	collector := []*sam.Record{}
	for a := range alignments {
		collector = append(collector, a)
	}
	if len(collector) != 2 {
		t.Fatalf("there should be 2 alignments")
	}
	for _, a := range collector {
		if a.Pos == 4 && a.Ref.Name() == "reference_A" {
			fmt.Println(a)
		} else if a.Pos == 3 && a.Ref.Name() == "reference_B" {
			fmt.Println(a)
		} else {
			t.Fatalf("alignment should be at position 4 of reference A or position 3 of reference B")
		}
	}
}

// this tests the second alignment method (seed shuffle)
func TestIncorrectSeed(t *testing.T) {
	// create the graph
	testGraph, testRefs, err := setupGraph(seq1, seq2)
	if err != nil {
		t.Fatal(err)
	}
	// create the read + incorrect seed
	testRead, err := setupRead(l3, l4, []byte("+"), []byte("===="), 2)
	if err != nil {
		t.Fatal(err)
	}
	// test the alignment
	alignments := Align(testRead, 0, &testGraph, testRefs, 5)
	collector := []*sam.Record{}
	for a := range alignments {
		collector = append(collector, a)
	}
	if len(collector) != 2 {
		t.Fatalf("there should be 2 alignments")
	}
	for _, a := range collector {
		if a.Pos == 4 && a.Ref.Name() == "reference_A" {
			fmt.Println(a)
		} else if a.Pos == 3 && a.Ref.Name() == "reference_B" {
			fmt.Println(a)
		} else {
			t.Fatalf("alignment should be at position 4 of reference A or position 3 of reference B")
		}
	}
}

// this tests the third alignment method (clipped read)
func TestClippedRead(t *testing.T) {
	// create the graph
	testGraph, testRefs, err := setupGraph(seq1, seq2)
	if err != nil {
		t.Fatal(err)
	}
	// create the read
	testRead, err := setupRead(l5, l6, []byte("+"), []byte("=========="), 0)
	if err != nil {
		t.Fatal(err)
	}
	// test the alignment
	alignments := Align(testRead, 0, &testGraph, testRefs, 5)
	for a := range alignments {
		fmt.Println(a)
		if a.Len() != 9 {
			t.Fatalf("alignment should have had 1 base clipped")
		}
	}
}
