package pipeline

import (
	"os"
	"testing"

	"github.com/will-rowe/groot/src/version"
)

///////////////////////////////////////////////////////////////////////////////////////////////

/*
TEST DATA
*/
// a single cluster of ~80 bla-OXA type genes
var msaList = []string{"test-data/test-genes.msa"}

// test reads derived from bla-OXA-90 and bla-OXA-106, both simulated with sequencing errors. OXA-106 has 4x as many reads as OXA-90
var fastq = []string{"test-data/test-reads-OXA90-OXA106-100bp-with-errors.fastq"}

// test reads derived from bla-OXA-90 only, simulated with sequencing errors
//var fastq = []string{"test-data/test-reads-OXA90-100bp-50x-with-errors.fastq"}

// the GFA produced by the sketch test
var gfaList = []string{"test-data/tmp/groot-graph-0.gfa"}

///////////////////////////////////////////////////////////////////////////////////////////////

/*
TEST PARAMETERS
*/
var testParameters = &Info{
	NumProc:              1,
	Version:              version.GetVersion(),
	KmerSize:             41,
	SketchSize:           50,
	WindowSize:           100,
	NumPart:              8,
	MaxK:                 4,
	MaxSketchSpan:        30,
	ContainmentThreshold: 0.99,
	IndexDir:             "test-data/tmp",
	Sketch: AlignCmd{
		MinKmerCoverage: 10,
		BloomFilter:     false,
		Fasta:           false,
		BAMout:          "test-data/tmp/out.bam",
	},
	Haplotype: HaploCmd{
		Cutoff:        1.0,
		MaxIterations: 10000,
		MinIterations: 50,
		HaploDir:      "test-data/tmp",
	},
}

func setupTmpDir() error {
	_ = os.RemoveAll("test-data/tmp")
	err := os.Mkdir("test-data/tmp", 0777)
	return err
}

///////////////////////////////////////////////////////////////////////////////////////////////

/*
DUMMY PIPELINE
*/

type ComponentA struct {
	input  []int
	output chan int
}

func NewComponentA(i []int) *ComponentA {
	return &ComponentA{input: i, output: make(chan int)}
}

func (ComponentA *ComponentA) Run() {
	defer close(ComponentA.output)
	for _, input := range ComponentA.input {
		ComponentA.output <- input
	}
}

type ComponentB struct {
	input    chan int
	addition int
	results  []int
}

func NewComponentB(i int) *ComponentB {
	return &ComponentB{addition: i}
}

func (ComponentB *ComponentB) Connect(previous *ComponentA) {
	ComponentB.input = previous.output
}

func (ComponentB *ComponentB) Run() {
	results := []int{}
	for input := range ComponentB.input {
		results = append(results, (input + ComponentB.addition))
	}
	ComponentB.results = results

}

///////////////////////////////////////////////////////////////////////////////////////////////

/*
DUMMY PIPELINE TEST
*/

func TestPipeline(t *testing.T) {
	inputValues := []int{1, 2, 3, 4}
	expectedOutput := []int{11, 12, 13, 14}

	// create the processes
	a := NewComponentA(inputValues)
	b := NewComponentB(10)

	// create the pipeline
	newPipeline := NewPipeline()

	// add the processes and connect them
	newPipeline.AddProcesses(a, b)
	b.Connect(a)
	if len(newPipeline.processes) != 2 {
		t.Fatal("did not add correct number of processes to pipeline")
	}

	// run the pipeline
	newPipeline.Run()

	// once the pipeline is done, there should be results in the final component
	if len(expectedOutput) != len(b.results) {
		t.Fatal("pipeline did not produce expected output")
	}
	for i, val := range b.results {
		if val != expectedOutput[i] {
			t.Fatal("pipeline did not produce expected output")
		}
	}
}
