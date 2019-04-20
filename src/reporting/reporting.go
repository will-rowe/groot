package reporting

import (
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

type annotation struct {
	arg      string
	count    int
	length   int
	coverage plotter.XYs
	cigar    string
}

type BAMreader struct {
	InputFile      string
	Plot           bool
	CoverageCutoff float64
	LowCov         bool
}

func NewBAMreader() *BAMreader {
	return &BAMreader{}
}

func (proc *BAMreader) Run() {
	// create a BAM reader from either STDIN or a BAM file
	var r io.Reader
	if proc.InputFile == "" {
		r = os.Stdin
	} else {
		f, err := os.Open(proc.InputFile)
		if err != nil {
			log.Fatalf("could not open BAM file %q:", err)
		}
		defer f.Close()
		ok, err := bgzf.HasEOF(f)
		if err != nil {
			log.Fatalf("could not open bam file %q:", err)
		}
		if !ok {
			log.Printf("file %q has no bgzf magic block: may be truncated", proc.InputFile)
		}
		r = f
	}
	b, err := bam.NewReader(r, 0)
	if err != nil {
		log.Fatalf("could not read BAM file: %q", err)
	}
	defer b.Close()

	// process the header
	argMap := make(map[string]*sam.Reference)
	for _, ref := range b.Header().Refs() {
		argMap[ref.Name()] = ref
	}

	// init a record map
	recordMap := make(map[string][]*sam.Record, len(argMap))
	for entry := range recordMap {
		recordMap[entry] = []*sam.Record{}
	}

	// process the records
	for {
		record, err := b.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		// TODO: this is where to skip records according to stringency of reporting (i.e. strict or moderate)
		// ignore unaligned
		if record.Flags == 4 {
			continue
		}
		// add the record to the corresponding reference sequence
		recordMap[record.Ref.Name()] = append(recordMap[record.Ref.Name()], record)
	}

	// launch a reporting goroutine for each reference sequence
	reportChan := make(chan annotation)
	var wg sync.WaitGroup
	for _, ref := range argMap {
		if records, ok := recordMap[ref.Name()]; ok {
			wg.Add(1)
			go func(recs []*sam.Record, ref *sam.Reference, sendChan chan<- annotation) {
				defer wg.Done()
				// coverageCheck tells us if all bases in the reference have been covered by a read
				coverageCheck := make(map[int]struct{})
				// pileup contains coverage value for each base in the reference
				pileup := make([]int, ref.Len())
				// for each record, move along the alignment and update reference coverage info
				for _, rec := range recs {
					recStart := rec.Start()
					recEnd := recStart + rec.Len()
					// if the read goes beyond the reference, only go up to the last base of the ref
					if recEnd > len(pileup)-1 {
						recEnd = len(pileup) - 1
					}

					for i := recStart; i <= recEnd; i++ {
						// if this is the first time the base has been covered, update coverageCheck
						if _, ok := coverageCheck[i]; !ok {
							coverageCheck[i] = struct{}{}
						}
						// update the pileup
						pileup[i]++
					}
				}

				// check we have a fully covered reference
				pileupCoverage := float64(len(coverageCheck)) / float64(len(pileup))
				if pileupCoverage >= proc.CoverageCutoff {
					// get the reference name (remove asterisk from cluster representative if it is present)
					refName := ref.Name()
					if refName[0] == 42 {
						refName = refName[1:]
					}
					// represent pileup as a CIGAR-ish string (so can see what bases aren't covered)
					cigar := []string{}
					// plot coverage for this gene using the pileup
					coverage := make(plotter.XYs, len(pileup))
					for i := range coverage {
						coverage[i].X = float64(i)
						coverage[i].Y = float64(pileup[i])
						if pileup[i] == 0 {
							cigar = append(cigar, "D")
						} else {
							cigar = append(cigar, "M")
						}
					}
					cleanCigar, internalD := cigarClean(cigar)
					if (internalD == true) && (proc.LowCov == true) {
						return
					}
					// create the annotation
					anno := annotation{
						arg:      refName,
						count:    len(records),
						length:   ref.Len(),
						coverage: coverage,
						cigar:    cleanCigar,
					}
					// send annotation on
					sendChan <- anno
				}
			}(records, ref, reportChan)
		}
	}
	go func() {
		wg.Wait()
		close(reportChan)
	}()

	// collect the annotated ARGs
	for anno := range reportChan {
		// print info to stdout
		fmt.Printf("%v\t%d\t%d\t%v\n", anno.arg, anno.count, anno.length, anno.cigar)

		// this will clean up the ARG name so that we can use it as a filename
		var replacer = strings.NewReplacer("/", "__", "\t", "__")

		// plot coverage for this gene
		if proc.Plot == true {
			covPlot, err := plot.New()
			if err != nil {
				panic(err)
			}
			covPlot.Title.Text = "coverage plot"
			covPlot.X.Label.Text = "position in gene"
			covPlot.Y.Label.Text = "coverage (number of reads at position)"
			err = plotutil.AddLinePoints(covPlot, anno.arg, anno.coverage)
			if err != nil {
				panic(err)
			}
			// clean the ARG name
			anno.arg = replacer.Replace(anno.arg)
			fileName := fmt.Sprintf("./groot-plots/coverage-for-%v.png", anno.arg)
			if err := covPlot.Save(8*vg.Inch, 8*vg.Inch, fileName); err != nil {
				panic(err)
			}
		}
	}
}

/*
  This function cleans up the cigar string
*/
func cigarClean(str []string) (string, bool) {
	counter := 1
	preVal := str[0]
	cigar := ""
	DMrecord := make(map[string]int)
	for i, val := range str {
		if i == 0 {
			continue
		}
		if i == len(str)-1 {
			if val == preVal {
				counter++
				cigar += strconv.Itoa(counter) + val
				DMrecord[val]++
			} else {
				cigar += strconv.Itoa(counter) + preVal + "1" + val
				DMrecord[val]++
			}
			break
		}
		if val == preVal {
			counter++
		} else {
			DMrecord[preVal]++
			cigar += strconv.Itoa(counter) + preVal
			preVal = val
			counter = 1
		}
	}
	// use the DM record to see if internal Ds have been found
	if ((DMrecord["D"] + DMrecord["M"]) <= 2) || ((DMrecord["D"] == 2) && (DMrecord["M"] == 1)) {
		return cigar, false
	} else {
		return cigar, true
	}
}
