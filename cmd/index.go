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
	"errors"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/pkg/profile"
	"github.com/spf13/cobra"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/lshForest"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/seqio"
	"github.com/will-rowe/groot/src/stream"
)

// the command line arguments
var (
	kSize        *int     // size of k-mer
	sigSize      *int     // size of MinHash signature
	readLength   *int     // length of query reads (used during alignment subcommand), needed as window length should ~= read length
	windowOffset *int     // used to move window through variation graph
	jsThresh     *float64 // minimum Jaccard similarity for LSH forest query
	msaDir       *string  // directory containing the input MSA files
	msaList      []string // the collected MSA files
	outDir       *string
)

// a default dir to store the index files
var defaultOutDir = "./groot-index-" + string(time.Now().Format("20060102150405"))

// the index command (used by cobra)
var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "Convert a set of clustered reference sequences to variation graphs and then index them",
	Long:  `Convert a set of clustered reference sequences to variation graphs and then index them`,
	Run: func(cmd *cobra.Command, args []string) {
		runIndex()
	},
}

/*
  A function to initialise the command line arguments
*/
func init() {
	RootCmd.AddCommand(indexCmd)
	kSize = indexCmd.Flags().IntP("kmerSize", "k", 7, "size of k-mer")
	sigSize = indexCmd.Flags().IntP("sigSize", "s", 128, "size of MinHash signature")
	readLength = indexCmd.Flags().IntP("readLength", "l", 100, "length of query reads (which will be aligned during the align subcommand)")
	windowOffset = indexCmd.Flags().IntP("windowOffset", "w", 1, "used to move window through variation graph")
	jsThresh = indexCmd.Flags().Float64P("jsThresh", "j", 0.99, "minimum Jaccard similarity for a seed to be recorded")
	msaDir = indexCmd.Flags().StringP("msaDir", "i", "", "directory containing the clustered references (MSA files)")
	outDir = indexCmd.PersistentFlags().StringP("outDir", "o", defaultOutDir, "directory to save index files to")
}

/*
  A function to check user supplied parameters
*/
func indexParamCheck() error {
	if *msaDir == "" {
		misc.ErrorCheck(errors.New("no MSA directory specified - run `groot index --help` for more info on the command"))
	}
	// check the we have received some MSA files TODO: could do with a better way of collecting these
	err := filepath.Walk(*msaDir, func(path string, f os.FileInfo, err error) error {
		if len(strings.Split(path, ".msa")) == 2 {
			msaList = append(msaList, path)
		}
		return nil
	})
	misc.ErrorCheck(err)
	if len(msaList) == 0 {
		return errors.New("no MSA files (.msa) found in the supplied directory")
	}
	// TODO: check the supplied arguments to make sure they don't conflict with each other eg:
	if *kSize > *readLength {
		return errors.New("supplied k-mer size greater than read length")
	}
	// setup the outDir
	if _, err := os.Stat(*outDir); os.IsNotExist(err) {
		if err := os.MkdirAll(*outDir, 0700); err != nil {
			return errors.New("can't create specified output directory")
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
	// set up logging
	logFH, err := os.OpenFile("groot-index.log", os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	defer logFH.Close()
	log.SetOutput(logFH)
	// set up profiling
	if *profiling == true {
		defer profile.Start(profile.ProfilePath("./")).Stop()
	}
	log.Printf("starting the index command")
	// check the supplied files and then log some stuff
	log.Printf("checking parameters...")
	misc.ErrorCheck(indexParamCheck())
	log.Printf("\tprocessors: %d", *proc)
	log.Printf("\tk-mer size: %d", *kSize)
	log.Printf("\tsignature size: %d", *sigSize)
	log.Printf("\tread length (window size): %d", *readLength)
	log.Printf("\twindow offset: %d", *windowOffset)
	log.Printf("\tnumber of MSA files found: %d", len(msaList))

	log.Printf("building graphs from MSAs...")
	// a store for the graphs
	graphStore := make(graph.GraphStore)
	skippedMSAs := 0
	// for each MSA, generate a variation graph
	for i, msa := range msaList {
		// create a pipeline
		pipeline := stream.NewPipeline()
		// Init processes
		dataStream := stream.NewDataStreamer()
		fastaHandler := stream.NewFastaHandler()
		msaChecker := stream.NewMsaChecker()
		// add in the user parameters
		dataStream.InputFile = []string{msa}
		msaChecker.MinSeqLength = *readLength
		// Arrange pipeline processes
		fastaHandler.Input = dataStream.Output
		msaChecker.Input = fastaHandler.Output
		// submit each process to the pipeline to be run
		pipeline.AddProcesses(dataStream, fastaHandler, msaChecker)
		// execute the pipeline in a go routine so that we can access the output channel of the final process (msaChecker)
		go func() {
			pipeline.Run()
		}()
		// store each sequence from the MSA in a slice
		args := []seqio.FASTAentry{}
		// if the MSA has failed the check, skip this msa
		discard := false
		for entry := range msaChecker.Output {
			if entry.Seq == nil {
				discard = true
			}
			args = append(args, entry)
		}
		// check the MSA was parsed
		if len(args) == 0 {
			misc.ErrorCheck(fmt.Errorf("could not parse \"%v\"", msa))
		}
		// check the MSA entries passed
		if discard == true {
			skippedMSAs++
			continue
		}
		// create nodes for the representative sequence
		graph, err := graph.NewGraph(args[0].Seq)
		misc.ErrorCheck(err)
		graph.Lengths = append(graph.Lengths, graph.NodeTotal)
		// add the ARG IDs to the graph
		for _, arg := range args {
			graph.ARGs = append(graph.ARGs, arg.ID)
		}
		// create nodes for the variation held in the other sequences
		if len(args) > 1 {
			for i := 1; i < len(args); i++ {
				before := graph.NodeTotal
				misc.ErrorCheck(graph.AddVariantNodes(args[i].Seq, i))
				// add the sequence length to the graph
				graph.Lengths = append(graph.Lengths, graph.NodeTotal-before)
			}
		}
		// store the graph
		graphStore[i-skippedMSAs] = graph
	} // finished parsing MSAs
	if len(graphStore) == 0 {
		misc.ErrorCheck(fmt.Errorf("could not create any graphs"))
	}
	log.Printf("\tnumber of graph built: %d", len(graphStore))
	log.Printf("windowing graphs and generating MinHash signatures...")
	// sort each graph and then generate signatures
	// create a channel to receive the signatures
	sigChan := make(chan graph.Window)
	// use go routines to sort graphs, window them and create signatures
	var wg sync.WaitGroup
	for ID, refGraph := range graphStore {
		wg.Add(1)
		go func(ID int, refGraph *graph.Graph) {
			defer wg.Done()
			// sort the graph
			misc.ErrorCheck(refGraph.TopSort())
			// iterate over each sequence held in the graph
			for i, j := 0, len(refGraph.ARGs); i < j; i++ {
				// get the length of the sequence
				windowingHalt := refGraph.Lengths[i] - *readLength
				// window the nodes belonging to the sequence in the graph
				for k := 0; k <= windowingHalt; k += *windowOffset {
					window, startNode, err := refGraph.Graph2Seq(i, k, *readLength)
					misc.ErrorCheck(err)
					// create a window and get MinHash signature
					newWindow := new(graph.Window)
					newWindow.Graph = ID
					newWindow.Node = startNode
					windowSeq := seqio.Sequence{Seq: window}
					newWindow.Sig = windowSeq.RunMinHash(*kSize, *sigSize).Signature()
					// send the window
					sigChan <- *newWindow
				}
			}
		}(ID, refGraph)
	}
	go func() {
		wg.Wait()
		close(sigChan)
	}()
	// create the sigStore to hold the signatures for each graph
	var sigStore = make([]map[int][][]uint64, len(graphStore))
	for i := range sigStore {
		sigStore[i] = make(map[int][][]uint64)
	}
	// receive the signatures and store
	var sigCount int = 0
	for window := range sigChan {
		sigStore[window.Graph][window.Node] = append(sigStore[window.Graph][window.Node], window.Sig)
		sigCount++
	}
	log.Printf("\tnumber of signatures generated: %d\n", sigCount)
	// run LSH forest
	log.Printf("running LSH forest...\n")
	database := lshForest.NewLSHforest(*sigSize, *jsThresh)
	// range over the nodes in each graph, each node will have one or more signature
	for graphID, nodesMap := range sigStore {
		// add each signature to the database
		for nodeID, signatures := range nodesMap {
			for _, signature := range signatures {
				// combine graphID and nodeID to form a string key for signature
				stringKey := fmt.Sprintf("g%dn%d", graphID, nodeID)
				// add the key to a lookup map
				key := seqio.Key{GraphID: graphID, Node: nodeID}
				database.KeyLookup[stringKey] = key
				// add the signature to the lshForest
				database.Add(stringKey, signature)
			}
		}
	}
	numHF, numBucks := database.Settings()
	log.Printf("\tnumber of hash functions per bucket: %d\n", numHF)
	log.Printf("\tnumber of buckets: %d\n", numBucks)
	// record runtime info
	info := &misc.IndexInfo{Ksize: *kSize, SigSize: *sigSize, JSthresh: *jsThresh, ReadLength: *readLength}
	// save the index files
	log.Printf("saving index files to \"%v\"...", *outDir)
	misc.ErrorCheck(info.Dump(*outDir + "/index.info"))
	log.Printf("\tsaved runtime info")
	misc.ErrorCheck(graphStore.Dump(*outDir + "/index.graph"))
	log.Printf("\tsaved variation graphs")
	misc.ErrorCheck(database.Dump(*outDir + "/index.sigs"))
	log.Printf("\tsaved MinHash signatures")
	log.Println("finished")
}
