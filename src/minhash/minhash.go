// Package minhash contains a minHash implementation that uses the ntHash rolling hash function
package minhash

import (
	"errors"
	"math"

	"github.com/will-rowe/ntHash"
)

// CANONICAL sets whether ntHash will return the canonical k-mer (inspects both strands for each k-mer and returns the lowest hash value)
const CANONICAL = false

// minHash struct contains all the minimum hash values for a sequence
type minHash struct {
	kSize     uint
	signature []uint64
}

// Add a sequence to the minHash
func (minHash *minHash) Add(sequence []byte) error {
	// initiate the rolling ntHash
	hasher, err := ntHash.New(&sequence, minHash.kSize)
	if err != nil {
		return err
	}
	// get hashed kmers from read
	var h1 uint32
	var h2 uint32
	for hv := range hasher.Hash(CANONICAL) {
		// for each hashed k-mer, try adding it to the sketch
		for i := 0; i < len(minHash.signature); i++ {
			h1 = uint32(hv)
			h2 = uint32(hv >> 32)
			newVal := uint64(h1) + (uint64(i) * uint64(h2))
			// evaluate and add to the signature if it is a minimum
			if newVal < minHash.signature[i] {
				minHash.signature[i] = newVal
			}
		}
	}
	return nil
}

// Signature returns the current minHash signature (set of minimums)
func (minHash *minHash) Signature() []uint64 {
	return minHash.signature
}

// Similarity is a method to estimate the Jaccard Similarity using to minHash signatures
func (minHash *minHash) Similarity(querySig []uint64) (float64, error) {
	if len(minHash.signature) != len(querySig) {
		return 0, errors.New("length of minhash signatures do not match")
	}
	intersect := 0
	for i := range minHash.signature {
		if minHash.signature[i] == querySig[i] {
			intersect++
		}
	}
	return float64(intersect) / float64(len(minHash.signature)), nil
}

// NewMinHash initiates a minHash struct and populates the signature with max values
func NewMinHash(kSize uint, sigSize int) *minHash {
	signature := make([]uint64, sigSize)
	for i := range signature {
		signature[i] = math.MaxUint64
	}
	return &minHash{
		kSize:     kSize,
		signature: signature,
	}
}
