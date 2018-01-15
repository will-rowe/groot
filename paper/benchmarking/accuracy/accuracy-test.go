package main

import (
	"flag"
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
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
		log.Printf("file %q has no bgzf magic block: may be truncated", inputFile)
	}
	r = f
	b, err := bam.NewReader(r, 0)
	if err != nil {
		log.Fatalf("could not read BAM file: %q", err)
	}
	defer b.Close()

	// process the records
	recordCount := 0
	readMap := make(map[string][]sam.Record)
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
		// populate readMap
		if _, ok := readMap[record.Name]; !ok {
			readMap[record.Name] = []sam.Record{*record}
		} else {
			readMap[record.Name] = append(readMap[record.Name], *record)
		}
		recordCount++
	}

	// iterate over reads, record false negatives and false positives
	unAligned, correctAligned, fPos, fNeg := *numTestReads, 0, 0, 0
	for read, hits := range readMap {
		unAligned--
		match := false
		// get the base range of the read
		splitReadName := strings.Split(read, "_argannot")
		splitReadName2 := strings.Split(splitReadName[0], "_")
		readPos := splitReadName2[len(splitReadName2)-1]
		// clean the read ID so it matches the reference ID
		readID := "argannot" + splitReadName[1]
		// make sure there is no whitespace in ID
		readID = strings.Split(readID, " ")[0]
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
				// see if the bases positions match
				if readPos != strconv.Itoa(hit.Pos) {
					fPos++
				} else {
					match = true
				}
			}
		}
		// record if this read has been correctly aligned
		if match == true {
			correctAligned++
		}
	}

	// output the numbers
	aligned := *numTestReads - unAligned
	al := float64(aligned) / float64(*numTestReads) * 100
	fmt.Printf("%v\t%.2f%%\t\taligned reads\n", aligned, al)
	ua := float64(unAligned) / float64(*numTestReads) * 100
	fmt.Printf("%v\t%.2f%%\t\tunaligned reads\n", unAligned, ua)
	misAligned := (*numTestReads - unAligned) - correctAligned
	ma := float64(misAligned) / float64(*numTestReads) * 100
	fmt.Printf("%d\t%.2f%%\t\tincorrectly aligned reads\n", misAligned, ma)
	// false negatives are any read that did not align or aligned wrongly (as the dataset used was generated from the database being queried)
	fNeg = unAligned + misAligned
	fn := float64(fNeg) / float64(*numTestReads) * 100
	fmt.Printf("%d\t%.2f%%\t\tfalse negatives (unaligned/incorrectly aligned READS)\n", fNeg, fn)
	// false positives refer to the number of incorrect ALIGNMENTS
	//fp := float64(fPos) / float64(recordCount) * 100
	//fmt.Printf("%d\t%.2f%%\t\tfalse positives (incorrect ALIGNMENTS)\n", fPos, fp)
}
