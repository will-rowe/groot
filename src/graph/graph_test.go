package graph

import (
	"io"
	"log"
	"os"
	"testing"

	"github.com/will-rowe/gfa"
)

var (
	inputFile  = "./test.gfa"
	inputFile2 = "./test.msa"
	windowSize = 150
	kmerSize   = 7
	sketchSize = 128
	blaB10     = []byte("ATGAAAGGATTAAAAGGGCTATTGGTTCTGGCTTTAGGCTTTACAGGACTACAGGTTTTTGGGCAACAGAACCCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAACCTATAATACCTTCAAAGGAACTAAATATGCGGCTAATGCGGTATATATGGTAACCGATAAAGGAGTAGTGGTTATAGACTCTCCATGGGGAGAAGATAAATTTAAAAGTTTTACAGACGAGATTTATAAAAAGCACGGAAAGAAAGTTATCATGAACATTGCAACCCACTCTCATGATGATAGAGCCGGAGGTCTTGAATATTTTGGTAAACTAGGTGCAAAAACTTATTCTACTAAAATGACAGATTCTATTTTAGCAAAAGAGAATAAGCCAAGAGCAAAGTACACTTTTGATAATAATAAATCTTTTAAAGTAGGAAAGACTGAGTTTCAGGTTTATTATCCGGGAAAAGGTCATACAGCAGATAATGTGGTTGTGTGGTTTCCTAAAGACAAAGTATTAGTAGGAGGCTGCATTGTAAAAAGTGGTGATTCGAAAGACCTTGGGTTTATTGGGGAAGCTTATGTAAACGACTGGACACAGTCCATACACAACATTCAGCAGAAATTTCCCTATGTTCAGTATGTCGTTGCAGGTCATGACGACTGGAAAGATCAAACATCAATACAACATACACTGGATTTAATCAGTGAATATCAACAAAAACAAAAGGCTTCAAATTAA")
)

func loadMSA() *gfa.GFA {
	// load the MSA
	msa, _ := gfa.ReadMSA(inputFile2)
	// convert the MSA to a GFA instance
	myGFA, err := gfa.MSA2GFA(msa)
	if err != nil {
		log.Fatal(err)
	}
	return myGFA
}

func loadGFA() *gfa.GFA {
	// load the GFA file
	fh, err := os.Open(inputFile)
	reader, err := gfa.NewReader(fh)
	if err != nil {
		log.Fatalf("can't read gfa file: %v", err)
	}
	// collect the GFA instance
	myGFA := reader.CollectGFA()
	// read the file
	for {
		line, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading line in gfa file: %v", err)
		}
		if err := line.Add(myGFA); err != nil {
			log.Fatalf("error adding line to GFA instance: %v", err)
		}
	}
	return myGFA
}

// test CreateGrootGraph
func TestCreateGrootGraph(t *testing.T) {
	myGFA, err := LoadGFA(inputFile)
	if err != nil {
		t.Fatal(err)
	}
	_, err = CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
}

// test Graph2Seq
func TestGraph2Seqs(t *testing.T) {
	t.Log("replace")
}

// test WindowGraph
func TestWindowGraph(t *testing.T) {
	myGFA := loadMSA()
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	counter := 0
	for window := range grootGraph.WindowGraph(windowSize, kmerSize, sketchSize) {
		//t.Log(window)
		_ = window
		counter++
	}
	t.Log("number of windows with unique sketchs: ", counter)
}

/*
// test ChainSegments
func TestChainSegments(t *testing.T) {
	myGFA := loadMSA()
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	chains, err := grootGraph.chainSegments(0, 3)
	if err != nil {
		t.Fatal(err)
	}
	counter := 0
	for _, chain := range chains {
		counter++
		t.Log(chain)
	}
	if counter != 2 {
		t.Fatal("two chains should be formed from this node")
	}
	err = grootGraph.BuildMarkovModel(1)
	if err != nil {
		t.Fatal(err)
	}
}
*/

// test FindMarkovPaths
func TestFindMarkovPaths(t *testing.T) {
	myGFA, err := LoadGFA("test2.gfa")
	if err != nil {
		t.Fatal(err)
	}
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	_ = grootGraph
}

// test SaveGraphAsGFA to save a gfa
func TestGraphDump(t *testing.T) {
	myGFA, err := LoadGFA(inputFile)
	if err != nil {
		t.Fatal(err)
	}
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	// add a dummy read so that the graph will write
	grootGraph.SortedNodes[0].IncrementKmerFreq(100.0)
	written, err := grootGraph.SaveGraphAsGFA("./tmp-graph.gfa", 0)
	if err != nil {
		t.Fatal(err)
	}
	if written != 1 {
		t.Fatal("graph not written as gfa file")
	}
	if err := os.Remove("./tmp-graph.gfa"); err != nil {
		t.Fatal(err)
	}
}
