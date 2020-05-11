// +build js,wasm

package bg

import (
	"fmt"
	"syscall/js"
	"time"

	"github.com/will-rowe/groot/src/lshe"
	"github.com/will-rowe/groot/src/pipeline"
)

// inputCheck is the callback to check the input is correct
func (GrootWASM *GrootWASM) inputCheck() interface{} {

	// check the input files first
	if len(GrootWASM.fastqFiles) == 0 {
		GrootWASM.statusUpdate("> no FASTQ files selected!")
		return nil
	}
	GrootWASM.iconUpdate("inputIcon")

	// check the parameters
	// TODO: this is where to check for alternative index requests
	GrootWASM.iconUpdate("paramIcon")

	// check the buffer has index data
	if len(GrootWASM.graphBuffer) == 0 {
		GrootWASM.statusUpdate("> no graphs in buffer")
		return nil
	}
	if len(GrootWASM.indexBuffer) == 0 {
		GrootWASM.statusUpdate("> no index in buffer")
		return nil
	}
	GrootWASM.inputChecked = true
	return nil
}

// setupGrootCb sets up the GROOT callback and runs GROOT when everything is set
func (GrootWASM *GrootWASM) setupGrootCb() {
	GrootWASM.grootCb = js.FuncOf(func(this js.Value, args []js.Value) interface{} {
		go func() {

			// stop GROOT?
			if GrootWASM.running == true {
				GrootWASM.running = false
				GrootWASM.statusUpdate("> stopped GROOT!")
				js.Global().Call("stopRecord")
				js.Global().Call("stopLogo")
				js.Global().Get("location").
					Call("reload")
				return
			}

			// check the input
			fmt.Println("checking input...")
			GrootWASM.inputCheck()
			if !GrootWASM.inputChecked {
				return
			}

			// start stuff
			js.Global().Call("startRecord")
			js.Global().Call("startLogo")
			GrootWASM.statusUpdate("> loading the index...")
			fmt.Println("loading index files...")

			// transfer the index from the buffer into a GROOT object
			if err := GrootWASM.info.LoadFromBytes(GrootWASM.graphBuffer); err != nil {
				GrootWASM.statusUpdate("> failed to load GROOT graphs!")
				fmt.Println(err)
				return
			}
			lshe := &lshe.ContainmentIndex{}
			if err := lshe.LoadFromBytes(GrootWASM.indexBuffer); err != nil {
				GrootWASM.statusUpdate("> failed to load GROOT index!")
				fmt.Println(err)
				return
			}
			GrootWASM.info.AttachDB(lshe)
			if GrootWASM.info == nil {
				return
			}

			// start GROOT notifications
			GrootWASM.running = true
			GrootWASM.statusUpdate("> running GROOT...")
			startTime := time.Now()

			// set up the pipeline
			sketchingPipeline := pipeline.NewPipeline()
			wasmStreamer := pipeline.NewWASMstreamer()
			fastqHandler := pipeline.NewFastqHandler(GrootWASM.info)
			fastqChecker := pipeline.NewFastqChecker(GrootWASM.info)
			readMapper := pipeline.NewReadMapper(GrootWASM.info)
			graphPruner := pipeline.NewGraphPruner(GrootWASM.info, true)
			emPathFinder := pipeline.NewEMpathFinder(GrootWASM.info)
			haploParser := pipeline.NewHaplotypeParser(GrootWASM.info)

			// connect the pipeline
			wasmStreamer.ConnectChan(GrootWASM.fastqInput)
			fastqHandler.ConnectWASM(wasmStreamer)
			fastqChecker.Connect(fastqHandler)
			readMapper.Connect(fastqChecker)
			graphPruner.Connect(readMapper)
			emPathFinder.ConnectPruner(graphPruner)
			haploParser.Connect(emPathFinder)
			sketchingPipeline.AddProcesses(wasmStreamer, fastqHandler, fastqChecker, readMapper, graphPruner, emPathFinder, haploParser)

			// start the stream and send data to the pipeline
			go js.Global().Call("fastqStreamer", GrootWASM.fastqFiles)

			// run the pipeline
			fmt.Println("starting the pipeline")
			sketchingPipeline.Run()
			fmt.Println("pipeline finished")

			// see how many reads mapped and number of k-mers processed
			readStats := readMapper.CollectReadStats()
			//foundPaths := graphPruner.CollectOutput()
			//foundHaplotypes := haploParser.CollectOutput()

			// print some updates
			if readStats[1] == 0 {
				fmt.Println("no reads mapped :(")
				js.Global().Call("stopRecord")
				js.Global().Call("stopLogo")
				GrootWASM.iconUpdate("startIcon")
				GrootWASM.statusUpdate("> no reads mapped to graphs :(")
				return
			}
			GrootWASM.statusUpdate(fmt.Sprintf("mapped reads = %d/%d", readStats[1], readStats[0]))

			// get the results
			GrootWASM.results = false
			for _, g := range GrootWASM.info.Store {
				paths, abundances := g.GetEMpaths()
				if len(paths) != 0 {
					GrootWASM.results = true
					fmt.Printf("\tgraph %d has %d called alleles after EM", g.GraphID, len(paths))
					for i, path := range paths {
						js.Global().Call("addResults", path, abundances[i])
					}
				}
			}

			// report any results
			js.Global().Call("stopRecord")
			js.Global().Call("stopLogo")
			GrootWASM.iconUpdate("startIcon")
			if GrootWASM.results == false {
				GrootWASM.statusUpdate("> no results found :(")
				fmt.Println("no paths left after running EM")
			} else {
				GrootWASM.statusUpdate("> GROOT finished!")
				secs := time.Since(startTime).Seconds()
				mins := time.Since(startTime).Minutes()
				timer := fmt.Sprintf("%.0fmins %.0fsecs", mins, secs)
				js.Global().Call("updateTimer", timer)
				js.Global().Call("showResults")
			}
		}()
		return nil
	})
}

//
func (GrootWASM *GrootWASM) printSeqs() {
	for _, g := range GrootWASM.info.Store {
		seqs, err := g.Graph2Seqs()
		if err != nil {
			GrootWASM.statusUpdate(fmt.Sprintf("%v", err))
		}
		for id, seq := range seqs {
			fmt.Printf(">%v\n%v\n", string(g.Paths[id]), string(seq))
		}
	}
}
