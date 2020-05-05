package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

var inputFile = flag.String("bam", "", "bam file to run accuracy test on")
var numTestReads = flag.Int("numReads", 0, "number of test reads used")

func main() {
	flag.Parse()
	var r io.Reader
	f, err := os.Open(*inputFile)
	if err != nil {
		log.Fatalf("could not open BAM file %q:", err)
	}
	defer f.Close()
	ok, err := bgzf.HasEOF(f)
	if err != nil {
		log.Fatalf("could not open bam file %q:", err)
	}
	if !ok {
		log.Printf("file %v has no bgzf magic block: may be truncated", inputFile)
	}
	r = f
	b, err := bam.NewReader(r, 0)
	if err != nil {
		log.Fatalf("could not read BAM file: %q", err)
	}
	defer b.Close()

	// process the records and keep all alignments in a map
	readMap := make(map[string][]sam.Record)
	multimapCount := 0
	for {
		record, err := b.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}

		// ignore unaligned
		if record.Flags == 4 {
			continue
		}

		// store alignments per read
		if _, ok := readMap[record.Name]; !ok {
			readMap[record.Name] = []sam.Record{*record}
		} else {
			if len(readMap[record.Name]) == 1 {
				multimapCount++
			}
			readMap[record.Name] = append(readMap[record.Name], *record)
		}
	}

	// work out the number of unaligned reads
	aligned := float64(len(readMap))
	unAligned := float64(*numTestReads) - aligned
	percAligned := aligned / float64(*numTestReads) * 100
	percMulti := float64(multimapCount) / float64(*numTestReads) * 100
	percUnaligned := float64(unAligned) / float64(*numTestReads) * 100
	fmt.Printf("%.0f\t%.2f%%\t\taligned reads\n", aligned, percAligned)
	fmt.Printf("%d\t%.2f%%\t\tmultialigned reads\n", multimapCount, percMulti)
	fmt.Printf("%.0f\t%.2f%%\t\tunaligned reads\n", unAligned, percUnaligned)

	// record false negatives and false positives
	correctAligned, correctStart, fPos := 0.0, 0.0, 0.0
	for read, hits := range readMap {
		match := false

		// get the read details from the read name (thanks to bbmap - randomreads.sh)
		splitReadName := strings.Split(read, "_")

		// get the read ID and remove whitespace so it matches the reference ID
		readID := strings.Split(splitReadName[9], "$")[0]
		readID = strings.Split(readID, " ")[0]

		// get the ref position the read was created from
		readRefPos, err := strconv.Atoi(splitReadName[2])
		if err != nil {
			log.Fatal(err)
		}

		// check all alignments for each read
		for _, hit := range hits {

			ref := hit.Ref.Name()

			// get rid of any leading * in the reference ID
			if hit.Ref.Name()[0] == 42 {
				ref = hit.Ref.Name()[1:]
			}

			// see if the alignment is to correct reference
			if ref != readID {
				fPos++
			} else {
				match = true

				// check that the aligment start pos is correct
				if hit.Start() != readRefPos {
					//fmt.Printf("start pos wrong for: %s (got %d, needed %d)\n", read, readRefPos, hit.Start())
				} else {
					correctStart++
				}
			}
		}

		// record if this read has been correctly aligned
		if match == true {
			correctAligned++
		}
	}

	// misaligned are reads that did not acheive their correct alignments
	misAligned := aligned - correctAligned
	percMisaligned := float64(misAligned) / float64(*numTestReads) * 100
	fmt.Printf("%.0f\t%.2f%%\t\tincorrectly aligned reads\n", misAligned, percMisaligned)

	// record the start pos stat - don't mark incorrect for now
	// fmt.Printf("\n%.0f/%.0f\tcorrect read alignments had correct start sites\n", correctStart, correctAligned)
}
