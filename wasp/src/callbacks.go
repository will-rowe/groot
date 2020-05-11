// +build js,wasm

package bg

import (
	"fmt"
	"strconv"
	"syscall/js"
)

// loadGraphs is linked with the JavaScript function of the same name
func (GrootWASM *GrootWASM) loadGraphs(this js.Value, args []js.Value) interface{} {
	fmt.Printf("loading graphs into buffer...\n")
	fmt.Printf("\tname: %v\n", args[0])
	fmt.Printf("\tsize: %d\n", args[2].Int())
	GrootWASM.graphBuffer = make([]byte, args[2].Int())
	size := js.CopyBytesToGo(GrootWASM.graphBuffer, args[1])
	if size != 0 {
		fmt.Printf("\tdone\n")
	} else {

		// shut down the app if any errors appeared
		GrootWASM.statusUpdate("> can't find graphs")
		GrootWASM.done <- struct{}{}
	}
	return nil
}

// loadIndex is linked with the JavaScript function of the same name
func (GrootWASM *GrootWASM) loadIndex(this js.Value, args []js.Value) interface{} {
	fmt.Printf("loading index into buffer...\n")
	fmt.Printf("\tname: %v\n", args[0])
	fmt.Printf("\tsize: %d\n", args[2].Int())
	GrootWASM.indexBuffer = make([]byte, args[2].Int())
	size := js.CopyBytesToGo(GrootWASM.indexBuffer, args[1])
	if size != 0 {
		fmt.Printf("\tdone\n")
	} else {

		// shut down the app if any errors appeared
		GrootWASM.statusUpdate("> can't find index")
		GrootWASM.done <- struct{}{}
	}
	return nil
}

// updateParameters is linked with the JavaScript function of the same name
func (GrootWASM *GrootWASM) updateParameters(this js.Value, args []js.Value) interface{} {
	fmt.Printf("updating parameters...\n")
	ct, err := strconv.ParseFloat(args[0].String(), 64)
	if err != nil {
		fmt.Println("error: could not parse float from parameter (cThresh)")
	}
	mk, err := strconv.ParseFloat(args[1].String(), 64)
	if err != nil {
		fmt.Println("error: could not parse float from parameter (mKmerCov)")
	}

	// shut down the app if any errors appeared
	if err != nil {
		GrootWASM.statusUpdate("> parameter parsing error")
		GrootWASM.done <- struct{}{}
		return nil
	}

	// otherwise, set all the parameters
	GrootWASM.info.ContainmentThreshold = ct
	GrootWASM.info.Sketch.MinKmerCoverage = mk

	// print a summary to the log
	fmt.Printf("\tcontainment theshold - %.2f\n", GrootWASM.info.ContainmentThreshold)
	fmt.Printf("\tminimum kmer coverage - %.0f\n", GrootWASM.info.Sketch.MinKmerCoverage)
	return nil
}

// getFiles is linked with the JavaScript function of the same name
//TODO: need to check this is correct
func (GrootWASM *GrootWASM) getFiles(this js.Value, args []js.Value) interface{} {
	fmt.Printf("loading fastq file list...\n")
	files := make([]interface{}, len(args))
	for i, val := range args {
		files[i] = val
	}
	if len(files) != 0 {
		fmt.Printf("\tadding files ready for WASM\n")
		GrootWASM.fastqFiles = append(GrootWASM.fastqFiles, files...)
	} else {
		fmt.Println("no input files found")
	}
	return nil
}

// munchFASTQ is linked with the JavaScript function of the same name
func (GrootWASM *GrootWASM) munchFASTQ(this js.Value, args []js.Value) interface{} {
	fastqBuffer := make([]byte, args[1].Int())
	_ = js.CopyBytesToGo(fastqBuffer, args[0])
	GrootWASM.fastqInput <- fastqBuffer
	return nil
}

// closeFASTQchan
func (GrootWASM *GrootWASM) closeFASTQchan(this js.Value, i []js.Value) interface{} {
	close(GrootWASM.fastqInput)
	return nil
}
