// Package pipeline contains a streaming pipeline implementation based on the Gopher Academy article by S. Lampa - Patterns for composable concurrent pipelines in Go (https://blog.gopheracademy.com/advent-2015/composable-pipelines-improvements/)
package pipeline

// BUFFERSIZE is the size of the buffer used by the pipeline channels
const BUFFERSIZE int = 64

// process is the interface used by pipeline
type process interface {
	Run()
}

// Pipeline is the base type, which takes any types that satisfy the process interface
type Pipeline struct {
	processes []process
}

// NewPipeline is the pipeline constructor
func NewPipeline() *Pipeline {
	return &Pipeline{}
}

// AddProcess is a method to add a single process to the pipeline
func (Pipeline *Pipeline) AddProcess(proc process) {
	// add the process to the pipeline
	Pipeline.processes = append(Pipeline.processes, proc)
}

// AddProcesses is a method to add multiple processes to the pipeline
func (Pipeline *Pipeline) AddProcesses(procs ...process) {
	for _, proc := range procs {
		Pipeline.AddProcess(proc)
	}
}

// Run is a method that starts the pipeline
func (Pipeline *Pipeline) Run() {
	// each pipeline process is run in a Go routines, except the last process which is run in the foreground to control the flow
	for i, process := range Pipeline.processes {
		if i < len(Pipeline.processes)-1 {
			go process.Run()
		} else {
			process.Run()
		}
	}
}

// GetNumProcesses is a method to return the number of processes registered in a pipeline
func (Pipeline *Pipeline) GetNumProcesses() int {
	return len(Pipeline.processes)
}
