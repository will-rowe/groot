package alignment

import (
	"io"
	"log"
	"os"
	"testing"

	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/seqio"
)

var (
	inputFile  = "../graph/test.gfa"
	windowSize = 100
	kSize      = 7
	sigSize    = 128
)

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

func setupMultimapRead() (*seqio.FASTQread, error) {
	testRead, err := seqio.NewFASTQread([]byte("@read-derived-from-segment-26"), []byte("CCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAAC"), []byte("+"), []byte("=================================================="))
	if err != nil {
		return nil, err
	}
	// create the seed
	seed := seqio.Key{
		GraphID: 1,
		Node:    26,
		OffSet:  0,
		RC:      false,
	}
	testRead.Seeds = append(testRead.Seeds, seed)
	return &testRead, nil
}

func setupUniqmapRead() (*seqio.FASTQread, error) {
	testRead, err := seqio.NewFASTQread([]byte("@read-derived-from-path-B10"), []byte("ATGAAAGGATTAAAAGGGCTATTGGTTCTGGCTTTAGGCTTTACAGGACTACAGGTTTTTGGGCAACAGAACCCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAACCTATAATACCTTCAAAGGAACTAAATATGCGGCTAATGCGGTATATATGGTAACCGATAAAGGAGTAGTGGTTATAGACTCTCCATGGGGAGAAGATAAATTTAAAAGTTTTACAGACGAGATTTATAAAAAGCACGGAAAGAAAGTTATCATGAACATTGCAACCCACTCTCATGATGATAGAGCCGGAGGTCTTGAATATTTTGGTAAACTAGGTGCAAAAACTTATTCTACTAAAATGACAGATTCTATTTTAGCAAAAGAGAATAAGCCAAGAGCAAAGTACACTTTTGATAATAATAAATCTTTTAAAGTAGGAAAGACTGAGTTTCAGGTTTATTATCCGGGAAAAGGTCATACAGCAGATAATGTGGTTGTGTGGTTTCCTAAAGACAAAGTATTAGTAGGAGGCTGCATTGTAAAAAGTGGTGATTCGAAAGACCTTGGGTTTATTGGGGAAGCTTATGTAAACGACTGGACACAGTCCATACACAACATTCAGCAGAAATTTCCCTATGTTCAGTATGTCGTTGCAGGTCATGACGACTGGAAAGATCAAACATCAATACAACATACACTGGATTTAATCAGTGAATATCAACAAAAACAAAAGGCTTCAAATTAA"), []byte("+"), []byte("ATGAAAGGATTAAAAGGGCTATTGGTTCTGGCTTTAGGCTTTACAGGACTACAGGTTTTTGGGCAACAGAACCCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAACCTATAATACCTTCAAAGGAACTAAATATGCGGCTAATGCGGTATATATGGTAACCGATAAAGGAGTAGTGGTTATAGACTCTCCATGGGGAGAAGATAAATTTAAAAGTTTTACAGACGAGATTTATAAAAAGCACGGAAAGAAAGTTATCATGAACATTGCAACCCACTCTCATGATGATAGAGCCGGAGGTCTTGAATATTTTGGTAAACTAGGTGCAAAAACTTATTCTACTAAAATGACAGATTCTATTTTAGCAAAAGAGAATAAGCCAAGAGCAAAGTACACTTTTGATAATAATAAATCTTTTAAAGTAGGAAAGACTGAGTTTCAGGTTTATTATCCGGGAAAAGGTCATACAGCAGATAATGTGGTTGTGTGGTTTCCTAAAGACAAAGTATTAGTAGGAGGCTGCATTGTAAAAAGTGGTGATTCGAAAGACCTTGGGTTTATTGGGGAAGCTTATGTAAACGACTGGACACAGTCCATACACAACATTCAGCAGAAATTTCCCTATGTTCAGTATGTCGTTGCAGGTCATGACGACTGGAAAGATCAAACATCAATACAACATACACTGGATTTAATCAGTGAATATCAACAAAAACAAAAGGCTTCAAATTAA"))
	if err != nil {
		return nil, err
	}
	// create the seed
	seed := seqio.Key{
		GraphID: 1,
		Node:    2,
		OffSet:  0,
		RC:      false,
	}
	testRead.Seeds = append(testRead.Seeds, seed)
	return &testRead, nil
}

// this tests the first alignment method (exact alignment) with a multimapper
func TestExactMatchMultiMapper(t *testing.T) {
	// create the read
	testRead, err := setupMultimapRead()
	if err != nil {
		log.Fatal(err)
	}
	// create the GrootGraph and graphStore
	myGFA := loadGFA()
	grootGraph, err := graph.CreateGrootGraph(myGFA, 1)
	if err != nil {
		log.Fatal(err)
	}
	graphStore := make(graph.GraphStore)
	graphStore[grootGraph.GraphID] = grootGraph
	referenceMap, err := graphStore.GetRefs()
	if err != nil {
		t.Fatal(err)
	}
	// align the read to the graph
	seedID := 0
	maxClip := 0
	alignments := Align(*testRead, seedID, grootGraph, referenceMap[grootGraph.GraphID], maxClip)
	count := 0
	for alignment := range alignments {
		t.Log(alignment)
		count++
	}
	if count != len(referenceMap[grootGraph.GraphID]) {
		t.Fatal("there should be an alignment returned for all references")
	}
}

// this tests the first alignment method (exact alignment) with a unique mapper
func TestExactMatchUniqMapper(t *testing.T) {
	// create the read
	testRead, err := setupUniqmapRead()
	if err != nil {
		log.Fatal(err)
	}
	// create the GrootGraph and graphStore
	myGFA := loadGFA()
	grootGraph, err := graph.CreateGrootGraph(myGFA, 1)
	if err != nil {
		log.Fatal(err)
	}
	graphStore := make(graph.GraphStore)
	graphStore[grootGraph.GraphID] = grootGraph
	referenceMap, err := graphStore.GetRefs()
	if err != nil {
		t.Fatal(err)
	}
	// align the read to the graph
	seedID := 0
	maxClip := 0
	alignments := Align(*testRead, seedID, grootGraph, referenceMap[grootGraph.GraphID], maxClip)
	count := 0
	for alignment := range alignments {
		t.Log(alignment)
		count++
	}
	if count != 1 {
		t.Fatal("there should be a single alignment")
	}
}
