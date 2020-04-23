/*
	the seqio package contains custom types and methods for holding and processing sequence data
*/
package seqio

import (
	"fmt"
	"unicode"

	"github.com/will-rowe/groot/src/minhash"
)

// encoding used by the FASTQ file
const encoding = 33

// complementBases is the lookup table used during reverse complementation
var complementBases = []byte{
	'A': 'T',
	'T': 'A',
	'C': 'G',
	'G': 'C',
	'N': 'N',
}

// Sequence is the base type for a FASTQ read
type Sequence struct {
	ID  []byte
	Seq []byte
}

// FASTQread is a type that holds a single FASTQ read, along with the locations it mapped to
type FASTQread struct {
	Sequence
	Misc []byte
	Qual []byte
	RC   bool
}

// RunMinHash is a method to create a minhash sketch for the sequence
func (Sequence *Sequence) RunMinHash(kmerSize, sketchSize int, kmv bool, bf *minhash.BloomFilter) ([]uint64, error) {

	// create the MinHash data structure, using the specified algorithm flavour
	var mh minhash.MinHash
	if kmv {
		mh = minhash.NewKMVsketch(uint(kmerSize), uint(sketchSize))
	} else {
		mh = minhash.NewKHFsketch(uint(kmerSize), uint(sketchSize))
	}

	// use the AddSequence method to populate the MinHash
	err := mh.AddSequence(Sequence.Seq)

	// get the sketch
	sketch := mh.GetSketch()

	// if the sketch isn't at capacity (in the case of BottomK sketches), fill up the remainder with 0s
	if kmv && len(sketch) != sketchSize {
		padding := make([]uint64, sketchSize-len(sketch))
		for i := 0; i < len(padding); i++ {
			padding[i] = 0
		}
		sketch = append(sketch, padding...)

	}

	// return the MinHash sketch and any error
	return sketch, err
}

// BaseCheck is a method to check for ACTGN bases and also to convert bases to upper case
// TODO: improve the efficiency of this...
func (Sequence *Sequence) BaseCheck() error {
	for i, j := 0, len(Sequence.Seq); i < j; i++ {
		switch base := unicode.ToUpper(rune(Sequence.Seq[i])); base {
		case 'A':
			Sequence.Seq[i] = byte(base)
		case 'C':
			Sequence.Seq[i] = byte(base)
		case 'T':
			Sequence.Seq[i] = byte(base)
		case 'G':
			Sequence.Seq[i] = byte(base)
		case 'N':
			Sequence.Seq[i] = byte(base)
		default:
			//return fmt.Errorf("non \"A\\C\\T\\G\\N\\-\" base (%v)", string(FASTQread.Seq[i]))
			Sequence.Seq[i] = byte('N')
		}
	}
	return nil
}

// RevComplement is a method to reverse complement a sequence held by a FASTQread
// TODO: the quality scores are currently not reversed by this method
func (FASTQread *FASTQread) RevComplement() {
	for i, j := 0, len(FASTQread.Seq); i < j; i++ {
		FASTQread.Seq[i] = complementBases[FASTQread.Seq[i]]
	}
	for i, j := 0, len(FASTQread.Seq)-1; i <= j; i, j = i+1, j-1 {
		FASTQread.Seq[i], FASTQread.Seq[j] = FASTQread.Seq[j], FASTQread.Seq[i]
	}
	if FASTQread.RC == true {
		FASTQread.RC = false
	} else {
		FASTQread.RC = true
	}
}

// QualTrim is a method to quality trim the sequence held by a FASTQread
/* the algorithm is based on bwa/cutadapt read quality trim functions:
-1. for each index position, subtract qual cutoff from the quality score
-2. sum these values across the read and trim at the index where the sum in minimal
-3. return the high-quality region
*/
func (FASTQread *FASTQread) QualTrim(minQual int) {
	start, qualSum, qualMax := 0, 0, 0
	end := len(FASTQread.Qual)
	for i, qual := range FASTQread.Qual {
		qualSum += minQual - (int(qual) - encoding)
		if qualSum < 0 {
			break
		}
		if qualSum > qualMax {
			qualMax = qualSum
			start = i + 1
		}
	}
	qualSum, qualMax = 0, 0
	for i, j := 0, len(FASTQread.Qual)-1; j >= i; j-- {
		qualSum += minQual - (int(FASTQread.Qual[j]) - encoding)
		if qualSum < 0 {
			break
		}
		if qualSum > qualMax {
			qualMax = qualSum
			end = j
		}
	}
	if start >= end {
		start, end = 0, 0
	}
	FASTQread.Seq = FASTQread.Seq[start:end]
	FASTQread.Qual = FASTQread.Qual[start:end]
}

// NewFASTQread generates a new fastq read from 4 lines of data
func NewFASTQread(l1 []byte, l2 []byte, l3 []byte, l4 []byte) (*FASTQread, error) {
	// check that it looks like a fastq read TODO: need more fastq checks
	//if len(l2) != len(l4) {
	//	return nil, errors.New("sequence and quality score lines are unequal lengths in fastq file")
	//}
	if l1[0] != 64 {
		return nil, fmt.Errorf("read ID in fastq file does not begin with @: %v", string(l1))
	}
	// create a FASTQread struct
	seq := Sequence{ID: l1, Seq: l2}
	return &FASTQread{
		Sequence: seq,
		Misc:     l3,
		Qual:     l4,
	}, nil
}
