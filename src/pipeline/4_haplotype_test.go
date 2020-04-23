package pipeline

import (
	"fmt"
	"os"
	"testing"
)

func TestHaplotyping(t *testing.T) {

	// load the files from the previous tests
	testParameters := new(Info)
	if err := testParameters.Load("test-data/tmp/groot.gg"); err != nil {
		t.Fatal(err)
	}
	haplotypingPipeline := NewPipeline()
	gfaReader := NewGFAreader(testParameters)
	emPathFinder := NewEMpathFinder(testParameters)
	haploParser := NewHaplotypeParser(testParameters)
	gfaReader.Connect(gfaList)
	emPathFinder.Connect(gfaReader)
	haploParser.Connect(emPathFinder)
	haplotypingPipeline.AddProcesses(gfaReader, emPathFinder, haploParser)
	if haplotypingPipeline.GetNumProcesses() != 3 {
		t.Fatal("wrong number of processes in pipeline")
	}
	haplotypingPipeline.Run()
	if len(testParameters.Store) != 1 {
		t.Fatal("haplotype should have recovered 1 graph")
	}

	for graphID, g := range testParameters.Store {
		fileName := fmt.Sprintf("test-data/tmp/groot-graph-%d-haplotype", graphID)
		_, err := g.SaveGraphAsGFA(fileName+".gfa", 0)
		if err != nil {
			t.Fatal(err)
		}
		seqs, err := g.Graph2Seqs()
		if err != nil {
			t.Fatal(err)
		}
		fh, err := os.Create(fileName + ".fna")
		if err != nil {
			t.Fatal(err)
		}
		for id, seq := range seqs {
			fmt.Fprintf(fh, ">%v\n%v\n", string(g.Paths[id]), string(seq))
		}
		fh.Close()
	}

	foundPaths := haploParser.CollectOutput()
	// this test needs improving - need to check for all the spiked in alleles
	correctPath := false
	for _, path := range foundPaths {
		t.Log(path)
		if path == "argannot~~~(Bla)OXA-90~~~EU547443:1-825" {
			correctPath = true
		}
	}
	if correctPath != true {
		t.Fatal("haplotyping did not identify correct allele in graph")
	}

	// remove the tmp files from all tests
	if err := os.Remove("test-data/tmp/groot.gg"); err != nil {
		t.Fatal("indexing did not create graph file: ", err)
	}
	if err := os.Remove("test-data/tmp/groot.lshe"); err != nil {
		t.Fatal("indexing did not create index file: ", err)
	}
	if err := os.Remove("test-data/tmp/groot-graph-0.gfa"); err != nil {
		t.Fatal("sketching did not create graph file: ", err)
	}
	if err := os.Remove("test-data/tmp/groot-graph-0-haplotype.fna"); err != nil {
		t.Fatal("haplotyping did not create fasta file: ", err)
	}
	if err := os.RemoveAll("test-data/tmp"); err != nil {
		t.Fatal("tests could not remove tmp directory")
	}

}
