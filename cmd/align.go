package cmd

import (
	"fmt"
	"log"
	"os"
	"runtime"
	"time"

	"github.com/pkg/profile"
	"github.com/spf13/cobra"
	"github.com/will-rowe/groot/src/lshe"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/pipeline"
	"github.com/will-rowe/groot/src/version"
)

// the command line arguments
var (
	fastq                *[]string                                                         // list of FASTQ files to align
	fasta                *bool                                                             // flag to treat input as fasta sequences
	containmentThreshold *float64                                                          // the containment threshold for the LSH ensemble
	minKmerCoverage      *float64                                                          // the minimum k-mer coverage per base of a segment
	graphDir             *string                                                           // directory to save gfa graphs to
	defaultGraphDir      = "./groot-graphs-" + string(time.Now().Format("20060102150405")) // a default graphDir
)

// alignCmd is used by cobra
var alignCmd = &cobra.Command{
	Use:   "align",
	Short: "Sketch sequences, align to references and weight variation graphs",
	Long:  `Sketch sequences, align to references and weight variation graphs`,
	Run: func(cmd *cobra.Command, args []string) {
		runSketch()
	},
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return misc.CheckRequiredFlags(cmd.Flags())
	},
}

// init the command line arguments
func init() {
	fastq = alignCmd.Flags().StringSliceP("fastq", "f", []string{}, "FASTQ file(s) to align")
	fasta = alignCmd.Flags().Bool("fasta", false, "if set, the input will be treated as fasta sequence(s) (experimental feature)")
	containmentThreshold = alignCmd.Flags().Float64P("contThresh", "t", 0.99, "containment threshold for the LSH ensemble")
	minKmerCoverage = alignCmd.Flags().Float64P("minKmerCov", "c", 1.0, "minimum number of k-mers covering each base of a graph segment")
	graphDir = alignCmd.PersistentFlags().StringP("graphDir", "g", defaultGraphDir, "directory to save variation graphs to")
	RootCmd.AddCommand(alignCmd)
}

// runAlign is the main function for the sketch sub-command
func runSketch() {

	// check index flag is set (global flag but don't require it for all sub commands)
	if *indexDir == "" {
		fmt.Println("please specify a directory with the index files (--indexDir)")
		os.Exit(1)
	}

	// set up profiling
	if *profiling {
		defer profile.Start(profile.MemProfile, profile.ProfilePath("./")).Stop()
		//defer profile.Start(profile.ProfilePath("./")).Stop()
	}

	// start logging
	if *logFile != "" {
		logFH := misc.StartLogging(*logFile)
		defer logFH.Close()
		log.SetOutput(logFH)
	} else {
		log.SetOutput(os.Stdout)
	}

	// start the sketch sub command
	start := time.Now()
	log.Printf("i am groot (version %s)", version.GetVersion())
	log.Printf("starting the sketch subcommand")

	// check the supplied files and then log some stuff
	log.Printf("checking parameters...")
	misc.ErrorCheck(alignParamCheck())
	log.Printf("\tminimum k-mer coverage: %.0f", *minKmerCoverage)
	log.Printf("\tprocessors: %d", *proc)
	for _, file := range *fastq {
		log.Printf("\tinput file: %v", file)
	}
	if *fasta {
		log.Print("\tinput file format: fasta")
	}
	log.Print("loading the index information...")
	info := new(pipeline.Info)
	misc.ErrorCheck(info.Load(*indexDir + "/groot.gg"))
	if info.Version != version.GetVersion() {
		misc.ErrorCheck(fmt.Errorf("the groot index was created with a different version of groot (you are currently using version %v)", version.GetVersion()))
	}
	log.Printf("\tk-mer size: %d\n", info.KmerSize)
	log.Printf("\tsketch size: %d\n", info.SketchSize)
	log.Printf("\twindow size used in indexing: %d\n", info.WindowSize)
	log.Print("loading the graphs...")
	log.Printf("\tnumber of variation graphs: %d\n", len(info.Store))
	log.Print("rebuilding the LSH Ensemble...")
	index := &lshe.ContainmentIndex{}
	misc.ErrorCheck(index.Load(*indexDir + "/groot.lshe"))
	info.AttachDB(index)
	if *profiling {
		log.Printf("\tloaded lshe file -> current memory usage %v", misc.PrintMemUsage())
		runtime.GC()
	}

	// add the sketch information to the existing groot runtime information
	info.NumProc = *proc
	info.Profiling = *profiling
	info.ContainmentThreshold = *containmentThreshold
	info.Sketch = pipeline.AlignCmd{
		Fasta:           *fasta,
		MinKmerCoverage: *minKmerCoverage,
	}
	log.Printf("\tcontainment threshold: %.2f\n", info.ContainmentThreshold)

	// create the pipeline
	log.Printf("initialising alignment pipeline...")
	alignmentPipeline := pipeline.NewPipeline()

	// initialise processes
	log.Printf("\tinitialising the processes")
	dataStream := pipeline.NewDataStreamer(info)
	fastqHandler := pipeline.NewFastqHandler(info)
	fastqChecker := pipeline.NewFastqChecker(info)
	readMapper := pipeline.NewReadMapper(info)
	graphPruner := pipeline.NewGraphPruner(info, false)

	// connect the pipeline processes
	log.Printf("\tconnecting data streams")
	dataStream.Connect(*fastq)
	fastqHandler.Connect(dataStream)
	fastqChecker.Connect(fastqHandler)
	readMapper.Connect(fastqChecker)
	graphPruner.Connect(readMapper)

	// submit each process to the pipeline and run it
	alignmentPipeline.AddProcesses(dataStream, fastqHandler, fastqChecker, readMapper, graphPruner)
	log.Printf("\tnumber of processes added to the alignment pipeline: %d\n", alignmentPipeline.GetNumProcesses())
	alignmentPipeline.Run()

	// once the sketching pipeline is finished, process the graph store and write the graphs to disk
	if len(info.Store) != 0 {
		log.Printf("saving graphs...\n")
		stats := readMapper.CollectReadStats()
		for graphID, g := range info.Store {
			fileName := fmt.Sprintf("%v/groot-graph-%d.gfa", *graphDir, graphID)
			_, err := g.SaveGraphAsGFA(fileName, stats[3])
			misc.ErrorCheck(err)
		}
	}
	log.Printf("finished in %s", time.Since(start))
}

// alignParamCheck is a function to check user supplied parameters
func alignParamCheck() error {

	// check the supplied FASTQ file(s)
	if len(*fastq) == 0 {
		misc.ErrorCheck(misc.CheckSTDIN())
		log.Printf("\tinput file: using STDIN")
	} else {
		for _, fastqFile := range *fastq {
			misc.ErrorCheck(misc.CheckFile(fastqFile))
			misc.ErrorCheck(misc.CheckExt(fastqFile, []string{"fastq", "fq", "fasta", "fna", "fa"}))
		}
	}

	// check the index directory and files
	misc.ErrorCheck(misc.CheckDir(*indexDir))
	misc.ErrorCheck(misc.CheckFile(*indexDir + "/groot.gg"))
	misc.ErrorCheck(misc.CheckFile(*indexDir + "/groot.lshe"))

	// setup the graphDir
	if _, err := os.Stat(*graphDir); os.IsNotExist(err) {
		if err := os.MkdirAll(*graphDir, 0700); err != nil {
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
