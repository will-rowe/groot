package graph

import (
	"fmt"
	"github.com/will-rowe/gfa"
	"io"
	"os"
	"testing"
)

var (
	inputFile  = "./test.gfa"
	inputFile2 = "./test.msa"
	windowSize = 150
	kSize      = 7
	sigSize    = 128
	blaB10     = []byte("ATGAAAGGATTAAAAGGGCTATTGGTTCTGGCTTTAGGCTTTACAGGACTACAGGTTTTTGGGCAACAGAACCCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAACCTATAATACCTTCAAAGGAACTAAATATGCGGCTAATGCGGTATATATGGTAACCGATAAAGGAGTAGTGGTTATAGACTCTCCATGGGGAGAAGATAAATTTAAAAGTTTTACAGACGAGATTTATAAAAAGCACGGAAAGAAAGTTATCATGAACATTGCAACCCACTCTCATGATGATAGAGCCGGAGGTCTTGAATATTTTGGTAAACTAGGTGCAAAAACTTATTCTACTAAAATGACAGATTCTATTTTAGCAAAAGAGAATAAGCCAAGAGCAAAGTACACTTTTGATAATAATAAATCTTTTAAAGTAGGAAAGACTGAGTTTCAGGTTTATTATCCGGGAAAAGGTCATACAGCAGATAATGTGGTTGTGTGGTTTCCTAAAGACAAAGTATTAGTAGGAGGCTGCATTGTAAAAAGTGGTGATTCGAAAGACCTTGGGTTTATTGGGGAAGCTTATGTAAACGACTGGACACAGTCCATACACAACATTCAGCAGAAATTTCCCTATGTTCAGTATGTCGTTGCAGGTCATGACGACTGGAAAGATCAAACATCAATACAACATACACTGGATTTAATCAGTGAATATCAACAAAAACAAAAGGCTTCAAATTAA")
)

func loadGFA() (*gfa.GFA, error) {
	// load the GFA file
	fh, err := os.Open(inputFile)
	reader, err := gfa.NewReader(fh)
	if err != nil {
		return nil, fmt.Errorf("can't read gfa file: %v", err)
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
			return nil, fmt.Errorf("error reading line in gfa file: %v", err)
		}
		if err := line.Add(myGFA); err != nil {
			return nil, fmt.Errorf("error adding line to GFA instance: %v", err)
		}
	}
	return myGFA, nil
}

func loadMSA() (*gfa.GFA, error) {
	// load the MSA
	msa, _ := gfa.ReadMSA(inputFile2)
	// convert the MSA to a GFA instance
	myGFA, err := gfa.MSA2GFA(msa)
	return myGFA, err
}

// test CreateGrootGraph
func TestCreateGrootGraph(t *testing.T) {
	myGFA, err := loadGFA()
	if err != nil {
		t.Fatal(err)
	}
	_, err = CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
}

// test Graph2Seq
func TestGraph2Seq(t *testing.T) {
	myGFA, err := loadGFA()
	if err != nil {
		t.Fatal(err)
	}
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	for pathID, pathName := range grootGraph.Paths {
		t.Log(string(pathName))
		t.Log(string(grootGraph.Graph2Seq(pathID)))
		if string(pathName) == "*argannot~~~(Bla)B-10~~~AY348325:1-747" {
			if string(grootGraph.Graph2Seq(pathID)) != string(blaB10) {
				t.Fatal("Graph2Seq did not reproduce blaB-10 sequence from GFA file")
			}
		}
	}
}

// test WindowGraph
func TestWindowGraph(t *testing.T) {
	myGFA, err := loadMSA()
	if err != nil {
		t.Fatal(err)
	}
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	counter := 0
	for window := range grootGraph.WindowGraph(windowSize, kSize, sigSize) {
		//t.Log(window)
		_ = window
		counter++
	}
	t.Log("number of windows with unique signatures: ", counter)
}

// test GraphStore dump/load
func TestGraphStore(t *testing.T) {
	myGFA, err := loadGFA()
	if err != nil {
		t.Fatal(err)
	}
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	graphStore := make(GraphStore)
	graphStore[0] = grootGraph
	if err := graphStore.Dump("./test.grootGraph"); err != nil {
		t.Fatal(err)
	}
	if err := graphStore.Load("./test.grootGraph"); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove("./test.grootGraph"); err != nil {
		t.Fatal(err)
	}
	// quick check of GetRefs
	if _, err := graphStore.GetRefs(); err != nil {
		t.Fatal(err)
	}
}

// test DumpGraph to save a gfa
func TestGraphDump(t *testing.T) {
	myGFA, err := loadGFA()
	if err != nil {
		t.Fatal(err)
	}
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		t.Fatal(err)
	}
	// add a dummy read so that the graph will write
	grootGraph.SortedNodes[0].IncrementReadCount()
	written, err := grootGraph.DumpGraph(".")
	if err != nil {
		t.Fatal(err)
	}
	if written != 1 {
		t.Fatal("graph not written as gfa file")
	}
	if err := os.Remove("*argannot~~~(Bla)B-10~~~AY3483251-747-groot-graph.gfa"); err != nil {
		t.Fatal(err)
	}
}
