package minhash

import (
	"fmt"
	"math"

	"github.com/will-rowe/nthash"
)

// KHFsketch is the structure for the K-Hash Functions MinHash sketch of a set of k-mers
type KHFsketch struct {
	kmerSize   uint
	sketchSize uint
	sketch     []uint64
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
func (KHFsketch *KHFsketch) AddSequence(sequence []byte) error {

	// initiate the rolling nthash
	hasher, err := nthash.NewHasher(&sequence, KHFsketch.kmerSize)
	if err != nil {
		return err
	}

	// range over the output of the hasher, where each iteration is a set of hash values for a k-mer
	for multiHashes := range hasher.MultiHash(CANONICAL, KHFsketch.sketchSize) {

		// evaluate if each hash value is lower than the existing one in the appropriate sketch position
		for i, min := range KHFsketch.sketch {
			if multiHashes[i] < min {
				KHFsketch.sketch[i] = multiHashes[i]
			}
		}
	}

	return nil
}

// GetSketch is a method to return the sketch held by a MinHash KHF sketch object
func (KHFsketch *KHFsketch) GetSketch() []uint64 {
	return KHFsketch.sketch
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
