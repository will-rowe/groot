package pipeline

import (
	"fmt"
	"os"
	"testing"

	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/misc"
)

func TestSketching(t *testing.T) {

	// load the files from the previous tests
	testParameters := new(Info)
	if err := testParameters.Load("test-data/tmp/groot.gg"); err != nil {
		t.Fatal(err)
	}

	lshe := &graph.ContainmentIndex{}
	if err := lshe.Load("test-data/tmp/groot.lshe"); err != nil {
		t.Fatal(err)
	}
	testParameters.AttachDB(lshe)

	// run the pipeline
	sketchingPipeline := NewPipeline()
	dataStream := NewDataStreamer(testParameters)
	fastqHandler := NewFastqHandler(testParameters)
	fastqChecker := NewFastqChecker(testParameters)
	readMapper := NewReadMapper(testParameters)
	graphPruner := NewGraphPruner(testParameters, false)
	dataStream.Connect(fastq)
	fastqHandler.Connect(dataStream)
	fastqChecker.Connect(fastqHandler)
	readMapper.Connect(fastqChecker)
	graphPruner.Connect(readMapper)
	sketchingPipeline.AddProcesses(dataStream, fastqHandler, fastqChecker, readMapper, graphPruner)
	if sketchingPipeline.GetNumProcesses() != 5 {
		t.Fatal("wrong number of processes in pipeline")
	}
	sketchingPipeline.Run()

	// check that the right number of reads mapped
	readStats := readMapper.CollectReadStats()
	t.Logf("total number of test reads = %d", readStats[0])
	t.Logf("number which mapped = %d", readStats[1])

	// check that we got the right allele in the approximately weighted graph
	foundPaths := graphPruner.CollectOutput()
	correctPath := false
	for _, path := range foundPaths {
		if path == "argannot~~~(Bla)OXA-90~~~EU547443:1-825" {
			correctPath = true
		}
	}
	if correctPath != true {
		t.Fatal("sketching did not identify correct allele in graph")
	}
	if err := testParameters.Dump("test-data/tmp/groot.gg"); err != nil {
		t.Fatal(err)
	}
	for graphID, g := range testParameters.Store {
		fileName := fmt.Sprintf("test-data/tmp/groot-graph-%d.gfa", graphID)
		_, err := g.SaveGraphAsGFA(fileName, readStats[3])
		misc.ErrorCheck(err)
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
	if err := os.RemoveAll("test-data/tmp"); err != nil {
		t.Fatal("tests could not remove tmp directory")
	}
}
