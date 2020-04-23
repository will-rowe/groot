package pipeline

import (
	"testing"
)

func TestIndexBuild(t *testing.T) {
	if err := setupTmpDir(); err != nil {
		t.Fatal(err)
	}
	indexingPipeline := NewPipeline()

	// initialise processes
	msaConverter := NewMSAconverter(testParameters)
	graphSketcher := NewGraphSketcher(testParameters)
	sketchIndexer := NewSketchIndexer(testParameters)

	// connect the pipeline processes
	msaConverter.Connect(msaList)
	graphSketcher.Connect(msaConverter)
	sketchIndexer.Connect(graphSketcher)

	// submit each process to the pipeline and run it
	indexingPipeline.AddProcesses(msaConverter, graphSketcher, sketchIndexer)
	if indexingPipeline.GetNumProcesses() != 3 {
		t.Fatal("wrong number of processes in pipeline")
	}
	indexingPipeline.Run()
	if err := testParameters.SaveDB("test-data/tmp/groot.lshe"); err != nil {
		t.Fatal(err)
	}
	if err := testParameters.Dump("test-data/tmp/groot.gg"); err != nil {
		t.Fatal(err)
	}
}

// benchmark indexing
func BenchmarkIndexing(b *testing.B) {
	// run the add method b.N times
	for n := 0; n < b.N; n++ {
		indexingPipeline := NewPipeline()
		msaConverter := NewMSAconverter(testParameters)
		graphSketcher := NewGraphSketcher(testParameters)
		sketchIndexer := NewSketchIndexer(testParameters)
		msaConverter.Connect(msaList)
		graphSketcher.Connect(msaConverter)
		sketchIndexer.Connect(graphSketcher)
		indexingPipeline.AddProcesses(msaConverter, graphSketcher, sketchIndexer)
		indexingPipeline.Run()
	}
}
