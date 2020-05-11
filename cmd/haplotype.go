package cmd

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"time"

	"github.com/pkg/profile"
	"github.com/spf13/cobra"
	"github.com/will-rowe/groot/src/graph"
	"github.com/will-rowe/groot/src/misc"
	"github.com/will-rowe/groot/src/pipeline"
	"github.com/will-rowe/groot/src/version"
)

// the command line arguments
var (
	graphDirectory  *string                                                              // directory containing the weighted variation graphs
	haploDir        *string                                                              // directory to write haplotype files to
	graphList       []string                                                             // the collected GFA files
	minIterations   *int                                                                 // minimum iterations for EM
	maxIterations   *int                                                                 // maximum iterations for EM
	cutOff          *float64                                                             // abundance cutoff for haplotypes
	probability     *float64                                                             // the probability threshold for a path to be reported
	defaultHaploDir = "./groot-haplotype-" + string(time.Now().Format("20060102150405")) // a default haploDir
)

// haplotypeCmd is used by cobra
var haplotypeCmd = &cobra.Command{
	Use:   "haplotype",
	Short: "Call haplotypes based on approximately weighted variation graphs",
	Long:  `Call haplotypes based on approximately weighted variation graphs`,
	Run: func(cmd *cobra.Command, args []string) {
		runHaplotype()
	},
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return misc.CheckRequiredFlags(cmd.Flags())
	},
}

// init the command line arguments
func init() {
	graphDirectory = haplotypeCmd.Flags().StringP("graphDirectory", "g", "", "directory containing the weighted variation graphs - required")
	haploDir = haplotypeCmd.PersistentFlags().StringP("outDir", "o", defaultHaploDir, "directory to write haplotype files to")
	minIterations = haplotypeCmd.PersistentFlags().IntP("minIterations", "x", 50, "minimum iterations for EM")
	maxIterations = haplotypeCmd.Flags().IntP("maxIterations", "y", 10000, "maximum iterations for EM")
	cutOff = haplotypeCmd.Flags().Float64P("cutOff", "z", 0.001, "abundance cutoff for calling haplotypes")
	haplotypeCmd.MarkFlagRequired("graphDirectory")
	RootCmd.AddCommand(haplotypeCmd)
}

// runHaplotype is the main function for the haplotype sub-command
func runHaplotype() {

	// check index flag is set (global flag but don't require it for all sub commands)
	if *indexDir == "" {
		fmt.Println("please specify a directory with the index files (--indexDir)")
		os.Exit(1)
	}

	// set up profiling
	if *profiling == true {
		//defer profile.Start(profile.MemProfile, profile.ProfilePath("./")).Stop()
		defer profile.Start(profile.ProfilePath("./")).Stop()
	}

	// start logging
	if *logFile != "" {
		logFH := misc.StartLogging(*logFile)
		defer logFH.Close()
		log.SetOutput(logFH)
	} else {
		log.SetOutput(os.Stdout)
	}
	// start sub command
	start := time.Now()
	log.Printf("i am groot (version %s)", version.GetVersion())
	log.Printf("starting the haplotype subcommand")
	// check the supplied files and then log some stuff
	log.Printf("checking parameters...")
	misc.ErrorCheck(haplotypeParamCheck())
	log.Printf("\tmin. iterations for EM: %d", *minIterations)
	log.Printf("\tmax. iterations for EM: %d", *maxIterations)
	log.Printf("\tabundance cut off reporting haplotypes: %0.2f", *cutOff)
	log.Printf("\tprocessors: %d", *proc)
	log.Print("loading the index...")
	log.Printf("\tgraph file: %v", *indexDir+"groot.gg")
	info := new(pipeline.Info)
	misc.ErrorCheck(info.Load(*indexDir + "/groot.gg"))
	if info.Version != version.GetVersion() {
		misc.ErrorCheck(fmt.Errorf("the groot index was created with a different version of groot (you are currently using version %v)", version.GetVersion()))
	}
	log.Printf("\tk-mer size: %d\n", info.KmerSize)
	log.Printf("\tsketch size: %d\n", info.SketchSize)
	log.Printf("\twindow size used in indexing: %d\n", info.WindowSize)
	log.Print("loading the graphs...")
	log.Printf("\tnumber of weighted GFAs for haplotyping: %d", len(graphList))

	// add the haplotype information to the existing groot runtime information
	info.Haplotype = pipeline.HaploCmd{
		MinIterations: *minIterations,
		MaxIterations: *maxIterations,
		Cutoff:        *cutOff,
		HaploDir:      *haploDir,
	}

	// create a graphStore
	info.Store = make(graph.Store)

	// create the pipeline
	log.Printf("initialising haplotype pipeline...")
	haplotypePipeline := pipeline.NewPipeline()

	// initialise processes
	log.Printf("\tinitialising the processes")
	gfaReader := pipeline.NewGFAreader(info)
	emPathFinder := pipeline.NewEMpathFinder(info)
	haploParser := pipeline.NewHaplotypeParser(info)

	// connect the pipeline processes
	log.Printf("\tconnecting data streams")
	gfaReader.Connect(graphList)
	emPathFinder.Connect(gfaReader)
	haploParser.Connect(emPathFinder)

	// submit each process to the pipeline and run it
	haplotypePipeline.AddProcesses(gfaReader, emPathFinder, haploParser)
	log.Printf("\tnumber of processes added to the haplotype pipeline: %d\n", haplotypePipeline.GetNumProcesses())
	log.Printf("finding alleles...")
	haplotypePipeline.Run()
	if len(info.Store) != 0 {
		log.Printf("writing files to \"%v/\"...\n", *haploDir)
		pathNames := []string{}
		for graphID, g := range info.Store {
			fileName := fmt.Sprintf("%v/groot-graph-%d-haplotype", *haploDir, graphID)
			_, err := g.SaveGraphAsGFA(fileName+".gfa", info.Haplotype.TotalKmers)
			misc.ErrorCheck(err)
			seqs, err := g.Graph2Seqs()
			misc.ErrorCheck(err)
			fh, err := os.Create(fileName + ".fna")
			misc.ErrorCheck(err)
			for id, seq := range seqs {
				fmt.Fprintf(fh, ">%v\n%v\n", string(g.Paths[id]), string(seq))
				pathNames = append(pathNames, string(g.Paths[id]))
			}
			fh.Close()
		}
		log.Printf("\tsaved haplotype graphs")
		log.Printf("\tsaved haplotype sequences")
		fileName := fmt.Sprintf("%v/summary.txt", *haploDir)
		fh, err := os.Create(fileName)
		misc.ErrorCheck(err)
		for _, name := range pathNames {
			fmt.Fprintf(fh, "%v\n", name)
		}
		fh.Close()
		log.Printf("\tsaved summary")
	}
	log.Printf("finished in %s", time.Since(start))
}

// haplotypeParamCheck is a function to check user supplied parameters
func haplotypeParamCheck() error {
	misc.ErrorCheck(misc.CheckDir(*indexDir))
	misc.ErrorCheck(misc.CheckFile(*indexDir + "/groot.gg"))
	misc.ErrorCheck(misc.CheckDir(*graphDirectory))
	graphs, err := filepath.Glob(*graphDirectory + "/groot-graph-*.gfa")
	if err != nil {
		return fmt.Errorf("can't find any graphs in the supplied graph directory")
	}
	for _, graph := range graphs {
		misc.ErrorCheck(misc.CheckFile(graph))
		graphList = append(graphList, graph)
	}
	// setup the haploDir
	if _, err := os.Stat(*haploDir); os.IsNotExist(err) {
		if err := os.MkdirAll(*haploDir, 0700); err != nil {
			return fmt.Errorf("can't create specified output directory")
		}
	}
	// check the cutoff
	if *cutOff > 1.0 || *cutOff < 0.0 {
		return fmt.Errorf("cutOff must be between 0.0 and 1.0")
	}
	// set number of processors to use
	if *proc <= 0 || *proc > runtime.NumCPU() {
		*proc = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*proc)
	return nil
}
