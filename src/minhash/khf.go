package minhash

import (
	"fmt"
	"math"
)

// HashValueSize is 8, the number of byte used for each hash value
const HashValueSize = 8
const seed = 42

// KHFsketch is the structure for the K-Hash Functions MinHash sketch of a set of k-mers
type KHFsketch struct {
	kmerSize   uint
	sketchSize uint
	sketch     []uint64
	hf1        func(b []byte) uint64
	hf2        func(b []byte) uint64
}

// NewKHFsketch is the constructor for a KHFsketch data structure
func NewKHFsketch(k, s uint) *KHFsketch {

	// init the sketch with maximum values
	sketch := make([]uint64, s)
	for i := range sketch {
		sketch[i] = math.MaxUint64
	}

	// return the data structure
	return &KHFsketch{
		kmerSize:   k,
		sketchSize: s,
		sketch:     sketch,
	}
}

// AddSequence is a method to decompose a read to canonical kmers, hash them and add any minimums to the sketch
func (mh *KHFsketch) AddSequence(sequence []byte) error {

	// check the sequence is long enough for given k
	if uint(len(sequence)) < mh.kmerSize {
		return fmt.Errorf("sequence length (%d) is short than k-mer length (%d)", len(sequence), mh.kmerSize)
	}

	// a holder for evalutating two k-mers
	kmers := [2]uint64{0, 0}

	// bitmask is used to update the previous k-mer with the next base
	bitmask := (uint64(1) << uint64(2*mh.kmerSize)) - uint64(1)
	bitshift := uint64(2 * (mh.kmerSize - 1))

	for i := 0; i < len(sequence); i++ {

		// get the nucleotide and convert to uint8
		c := seqNT4table[sequence[i]]

		// if the nucleotide == N
		if c > 3 {

			// TODO: handle these

		}

		// get the forward k-mer
		kmers[0] = (kmers[0]<<2 | uint64(c)) & bitmask

		// get the reverse k-mer
		kmers[1] = (kmers[1] >> 2) | (uint64(3)-uint64(c))<<bitshift

		// get the span of the k-mer
		if uint(i+1) < mh.kmerSize {
			continue
		}

		// set the canonical k-mer
		var strand uint
		if kmers[0] > kmers[1] {
			strand = 1
		}

		// get the two base hashes
		//hv1 := hash64(kmers[strand], bitmask)<<8 | uint64(mh.kmerSize)
		//hv2 := splitmix64(kmers[strand])
		hv2 := hash64(kmers[strand], bitmask)<<8 | uint64(mh.kmerSize)

		// try adding the k-mer in each slot of the sketch
		for j, min := range mh.sketch {
			hv := kmers[strand] + uint64(j)*hv2
			if hv < min {
				mh.sketch[j] = hv
			}
		}
	}

	return nil
}

// GetSketch is a method to return the sketch held by a MinHash KHF sketch object
func (mh *KHFsketch) GetSketch() []uint64 {
	return mh.sketch
}

// GetSimilarity estimates the similarity between two k-mer sets based on the KHF sketch
func (mh *KHFsketch) GetSimilarity(mh2 MinHash) (float64, error) {

	// check this is a pair of KHF
	if fmt.Sprintf("%T", mh) != fmt.Sprintf("%T", mh2) {
		return 0.0, fmt.Errorf("mismatched MinHash types: %T vs. %T", mh, mh2)
	}
	querySketch, ok := mh2.(*KHFsketch)
	if !ok {
		return 0.0, fmt.Errorf("could not assert sketch is a KHF")
	}

	// check sketch lengths match
	if len(mh.sketch) != len(querySketch.sketch) {
		return 0.0, fmt.Errorf("sketches do not have the same number of minimums: %d vs %d", len(mh.sketch), len(querySketch.sketch))
	}

	// find intersection
	intersect := 0
	for i := uint(0); i < mh.sketchSize; i++ {
		if mh.sketch[i] == querySketch.sketch[i] {
			intersect++
		}
	}

	// return jaccard similarity estimate
	return float64(intersect) / float64(mh.sketchSize), nil
}
