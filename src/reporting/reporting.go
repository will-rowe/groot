package reporting

import (
	"fmt"
	"io"
	"log"
	"os"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

type annotation struct {
	arg    string
	count  int
	length int
}

// TODO - add total number of reads - either do this with a SAM header field or will have to include unmapped reads in SAM

type BAMreader struct {
	InputFile      string
	CoverageCutoff float64
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
					// get the reference name (remove asterisk if present - from cluster representative)
					refName := ref.Name()
					if refName[0] == 42 {
						refName = refName[1:]
					}
					// create the annotation
					anno := annotation{
						arg:    refName,
						count:  len(records),
						length: ref.Len(),
					}
					// send it on
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
		fmt.Printf("%v\t%d\t%d\n", anno.arg, anno.count, anno.length)
	}
}
