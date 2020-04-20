// Copyright © 2017 Will Rowe <will.rowe@stfc.ac.uk>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package cmd

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/biogo/biogo/seq/multi"
	"github.com/pkg/profile"
	"github.com/spf13/cobra"
	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/lshForest"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
	"github.com/will-rowe/groot/src/version"
)

// the command line arguments
var (
	kSize         *uint                                                            // size of k-mer
	sigSize       *int                                                             // size of MinHash signature
	readLength    *int                                                             // length of query reads (used during alignment subcommand), needed as window length should ~= read length
	jsThresh      *float64                                                         // minimum Jaccard similarity for LSH forest query
	msaDir        *string                                                          // directory containing the input MSA files
	msaList       []string                                                         // the collected MSA files
	outDir        *string                                                          // directory to save index files and log to
	defaultOutDir = "./groot-index-" + string(time.Now().Format("20060102150405")) // a default dir to store the index files
)

// the index command (used by cobra)
var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "Convert a set of clustered reference sequences to variation graphs and then index them",
	Long:  `Convert a set of clustered reference sequences to variation graphs and then index them`,
	Run: func(cmd *cobra.Command, args []string) {
		runIndex()
	},
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return misc.CheckRequiredFlags(cmd.Flags())
	},
}

// a function to initialise the command line arguments
func init() {
	kSize = indexCmd.Flags().UintP("kmerSize", "k", 7, "size of k-mer")
	sigSize = indexCmd.Flags().IntP("sigSize", "s", 128, "size of MinHash signature")
	readLength = indexCmd.Flags().IntP("readLength", "l", 100, "length of query reads (which will be aligned during the align subcommand)")
	jsThresh = indexCmd.Flags().Float64P("jsThresh", "j", 0.99, "minimum Jaccard similarity for a seed to be recorded")
	msaDir = indexCmd.Flags().StringP("msaDir", "i", "", "directory containing the clustered references (MSA files) - required")
	outDir = indexCmd.PersistentFlags().StringP("outDir", "o", defaultOutDir, "directory to save index files to")
	indexCmd.MarkFlagRequired("msaDir")
	RootCmd.AddCommand(indexCmd)
}

//  a function to check user supplied parameters
func indexParamCheck() error {
	if *msaDir == "" {
		misc.ErrorCheck(fmt.Errorf("no MSA directory specified - run `groot index --help` for more info on the command"))
	}
	if _, err := os.Stat(*msaDir); os.IsNotExist(err) {
		return fmt.Errorf("can't find specified MSA directory")
	}
	// check the we have received some MSA files TODO: could do with a better way of collecting these
	err := filepath.Walk(*msaDir, func(path string, f os.FileInfo, err error) error {
		// ignore dot files
		if f.Name()[0] == 46 {
			return nil
		}
		// ignore empty files
		if f.Size() == 0 {
			return nil
		}
		// keep anything with an .msa extension
		if len(strings.Split(path, ".msa")) == 2 {
			msaList = append(msaList, path)
		}
		return nil
	})
	misc.ErrorCheck(err)
	if len(msaList) == 0 {
		return fmt.Errorf("no MSA files (.msa) found in the supplied directory")
	}
	// TODO: check the supplied arguments to make sure they don't conflict with each other eg:
	if *kSize > uint(*readLength) {
		return fmt.Errorf("supplied k-mer size greater than read length")
	}
	// setup the outDir
	if _, err := os.Stat(*outDir); os.IsNotExist(err) {
		if err := os.MkdirAll(*outDir, 0700); err != nil {
			return fmt.Errorf("can't create specified output directory")
		}
	}
	// set number of processors to use
	if *proc <= 0 || *proc > runtime.NumCPU() {
		*proc = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*proc)
	return nil
}

/*
  The main function for the index command
*/
func runIndex() {
	// set up profiling
	if *profiling == true {
		defer profile.Start(profile.ProfilePath("./")).Stop()
	}
	// start logging
	logFH := misc.StartLogging(*logFile)
	defer logFH.Close()
	log.SetOutput(logFH)
	log.Printf("i am groot (version %s)", version.VERSION)
	log.Printf("starting the index subcommand")
	// check the supplied files and then log some stuff
	log.Printf("checking parameters...")
	misc.ErrorCheck(indexParamCheck())
	log.Printf("\tprocessors: %d", *proc)
	log.Printf("\tk-mer size: %d", *kSize)
	log.Printf("\tsignature size: %d", *sigSize)
	log.Printf("\tread length (window size): %d", *readLength)
	log.Printf("\tnumber of MSA files found: %d", len(msaList))
	///////////////////////////////////////////////////////////////////////////////////////
	log.Printf("building groot graphs...")
	// process each msa in a go routine
	var wg sync.WaitGroup
	graphChan := make(chan *graph.GrootGraph)
	for i, msaFile := range msaList {
		// load the MSA outside of the go-routine to prevent 'too many open files' error on OSX
		msa, err := gfa.ReadMSA(msaFile)
		misc.ErrorCheck(err)
		wg.Add(1)
		go func(msaID int, msa *multi.Multi) {
			defer wg.Done()
			// convert the MSA to a GFA instance
			newGFA, err := gfa.MSA2GFA(msa)
			misc.ErrorCheck(err)
			// create a GrootGraph
			grootGraph, err := graph.CreateGrootGraph(newGFA, msaID)
			if err != nil {
				log.Fatal(err)
			}
			graphChan <- grootGraph
		}(i, msa)
	}
	go func() {
		wg.Wait()
		close(graphChan)
	}()
	///////////////////////////////////////////////////////////////////////////////////////
	// collect and store the GrootGraphs
	graphStore := make(graph.GraphStore)
	for graph := range graphChan {
		graphStore[graph.GraphID] = graph
	}
	if len(graphStore) == 0 {
		misc.ErrorCheck(fmt.Errorf("could not create any graphs"))
	}
	log.Printf("\tnumber of groot graphs built: %d", len(graphStore))
	///////////////////////////////////////////////////////////////////////////////////////
	log.Printf("windowing graphs and generating MinHash signatures...")
	// process each graph in a go routine
	windowChan := make(chan *graph.Window)
	for _, grootGraph := range graphStore {
		wg.Add(1)
		go func(grootGraph *graph.GrootGraph) {
			defer wg.Done()
			// create signature for each window in the graph
			for window := range grootGraph.WindowGraph(*readLength, *kSize, *sigSize) {
				windowChan <- window
			}
		}(grootGraph)
	}
	go func() {
		wg.Wait()
		close(windowChan)
	}()
	///////////////////////////////////////////////////////////////////////////////////////
	// collect and store the GrootGraph windows
	var sigStore = make([]map[int]map[int][][]uint64, len(graphStore))
	for i := range sigStore {
		sigStore[i] = make(map[int]map[int][][]uint64)
	}
	// receive the signatures
	var sigCount int = 0
	for window := range windowChan {
		// initialise the inner map of sigStore if graph has not been seen yet
		if _, ok := sigStore[window.GraphID][window.Node]; !ok {
			sigStore[window.GraphID][window.Node] = make(map[int][][]uint64)
		}
		// store the signature for the graph:node:offset
		sigStore[window.GraphID][window.Node][window.OffSet] = append(sigStore[window.GraphID][window.Node][window.OffSet], window.Sig)
		sigCount++
	}
	log.Printf("\tnumber of signatures generated: %d\n", sigCount)
	///////////////////////////////////////////////////////////////////////////////////////
	// run LSH forest
	log.Printf("running LSH forest...\n")
	database := lshForest.NewLSHforest(*sigSize, *jsThresh)
	// range over the nodes in each graph, each node will have one or more signature
	for graphID, nodesMap := range sigStore {
		// add each signature to the database
		for nodeID, offsetMap := range nodesMap {
			for offset, signatures := range offsetMap {
				for _, signature := range signatures {
					// combine graphID, nodeID and offset to form a string key for signature
					stringKey := fmt.Sprintf("g%dn%do%d", graphID, nodeID, offset)
					// add the key to a lookup map
					key := seqio.Key{GraphID: graphID, Node: nodeID, OffSet: offset}
					database.KeyLookup[stringKey] = key
					// add the signature to the lshForest
					database.Add(stringKey, signature)
				}
			}
		}
	}
	numHF, numBucks := database.Settings()
	log.Printf("\tnumber of hash functions per bucket: %d\n", numHF)
	log.Printf("\tnumber of buckets: %d\n", numBucks)
	///////////////////////////////////////////////////////////////////////////////////////
	// record runtime info
	info := &misc.IndexInfo{Version: version.VERSION, Ksize: *kSize, SigSize: *sigSize, JSthresh: *jsThresh, ReadLength: *readLength}
	// save the index files
	log.Printf("saving index files to \"%v\"...", *outDir)
	misc.ErrorCheck(info.Dump(*outDir + "/index.info"))
	log.Printf("\tsaved runtime info")
	misc.ErrorCheck(graphStore.Dump(*outDir + "/index.graph"))
	log.Printf("\tsaved groot graphs")
	misc.ErrorCheck(database.Dump(*outDir + "/index.sigs"))
	log.Printf("\tsaved MinHash signatures")
	log.Println("finished")
}
