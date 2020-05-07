package graph

import (
	"log"
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/groot/src/lshe"
	"github.com/will-rowe/groot/src/seqio"
)

func setupRead() (*seqio.FASTQread, *lshe.Key, error) {
	testRead, err := seqio.NewFASTQread([]byte("@read-derived-from-path-B7"), []byte("ATGAAAGGATTAAAAGGG"), []byte("+"), []byte("++++++++++++++++++"))
	if err != nil {
		return nil, nil, err
	}
	seed := &lshe.Key{
		GraphID: 1,
		Node:    2,
		OffSet:  0,
		RC:      false,
	}
	return testRead, seed, nil
}

func setupMultimapRead() (*seqio.FASTQread, *lshe.Key, error) {
	testRead, err := seqio.NewFASTQread([]byte("@read-derived-from-segment-26"), []byte("CCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAAC"), []byte("+"), []byte("=================================================="))
	if err != nil {
		return nil, nil, err
	}
	seed := &lshe.Key{
		GraphID: 1,
		Node:    26,
		OffSet:  0,
		RC:      false,
	}
	return testRead, seed, nil
}

func setupUniqmapRead() (*seqio.FASTQread, *lshe.Key, error) {
	testRead, err := seqio.NewFASTQread([]byte("@read-derived-from-path-B10"), []byte("ATGAAAGGATTAAAAGGGCTATTGGTTCTGGCTTTAGGCTTTACAGGACTACAGGTTTTTGGGCAACAGAACCCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAACCTATAATACCTTCAAAGGAACTAAATATGCGGCTAATGCGGTATATATGGTAACCGATAAAGGAGTAGTGGTTATAGACTCTCCATGGGGAGAAGATAAATTTAAAAGTTTTACAGACGAGATTTATAAAAAGCACGGAAAGAAAGTTATCATGAACATTGCAACCCACTCTCATGATGATAGAGCCGGAGGTCTTGAATATTTTGGTAAACTAGGTGCAAAAACTTATTCTACTAAAATGACAGATTCTATTTTAGCAAAAGAGAATAAGCCAAGAGCAAAGTACACTTTTGATAATAATAAATCTTTTAAAGTAGGAAAGACTGAGTTTCAGGTTTATTATCCGGGAAAAGGTCATACAGCAGATAATGTGGTTGTGTGGTTTCCTAAAGACAAAGTATTAGTAGGAGGCTGCATTGTAAAAAGTGGTGATTCGAAAGACCTTGGGTTTATTGGGGAAGCTTATGTAAACGACTGGACACAGTCCATACACAACATTCAGCAGAAATTTCCCTATGTTCAGTATGTCGTTGCAGGTCATGACGACTGGAAAGATCAAACATCAATACAACATACACTGGATTTAATCAGTGAATATCAACAAAAACAAAAGGCTTCAAATTAA"), []byte("+"), []byte("ATGAAAGGATTAAAAGGGCTATTGGTTCTGGCTTTAGGCTTTACAGGACTACAGGTTTTTGGGCAACAGAACCCTGATATTAAAATTGAAAAATTAAAAGATAATTTATACGTCTATACAACCTATAATACCTTCAAAGGAACTAAATATGCGGCTAATGCGGTATATATGGTAACCGATAAAGGAGTAGTGGTTATAGACTCTCCATGGGGAGAAGATAAATTTAAAAGTTTTACAGACGAGATTTATAAAAAGCACGGAAAGAAAGTTATCATGAACATTGCAACCCACTCTCATGATGATAGAGCCGGAGGTCTTGAATATTTTGGTAAACTAGGTGCAAAAACTTATTCTACTAAAATGACAGATTCTATTTTAGCAAAAGAGAATAAGCCAAGAGCAAAGTACACTTTTGATAATAATAAATCTTTTAAAGTAGGAAAGACTGAGTTTCAGGTTTATTATCCGGGAAAAGGTCATACAGCAGATAATGTGGTTGTGTGGTTTCCTAAAGACAAAGTATTAGTAGGAGGCTGCATTGTAAAAAGTGGTGATTCGAAAGACCTTGGGTTTATTGGGGAAGCTTATGTAAACGACTGGACACAGTCCATACACAACATTCAGCAGAAATTTCCCTATGTTCAGTATGTCGTTGCAGGTCATGACGACTGGAAAGATCAAACATCAATACAACATACACTGGATTTAATCAGTGAATATCAACAAAAACAAAAGGCTTCAAATTAA"))
	if err != nil {
		return nil, nil, err
	}
	seed := &lshe.Key{
		GraphID: 1,
		Node:    2,
		OffSet:  0,
		RC:      false,
	}
	return testRead, seed, nil
}

func setupGraph() (*GrootGraph, map[int][]*sam.Reference, error) {
	myGFA := loadGFA()
	grootGraph, err := CreateGrootGraph(myGFA, 1)
	if err != nil {
		return nil, nil, err
	}
	graphStore := make(Store)
	graphStore[grootGraph.GraphID] = grootGraph
	references, err := graphStore.GetSAMrefs()
	if err != nil {
		return nil, nil, err
	}
	return grootGraph, references, nil
}

// this tests the first alignment method (exact alignment) with a multimapper
func TestExactMatchMultiMapper(t *testing.T) {

	// get the read and multimap seed
	testRead, seed, err := setupRead()
	if err != nil {
		log.Fatal(err)
	}

	// create the GrootGraph and SAM refs for this graph
	grootGraph, references, err := setupGraph()
	if err != nil {
		t.Fatal(err)
	}

	// align the read to the graph
	alignments, err := grootGraph.AlignRead(testRead, seed, references[int(grootGraph.GraphID)])
	if err != nil {
		t.Fatal(err)
	}
	for _, alignment := range alignments {
		t.Log(alignment.String())
		t.Log(alignment.Start())
	}

}
