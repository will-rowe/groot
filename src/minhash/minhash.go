// the minhash package contains a minHash implementation that uses the ntHash rolling hash function
package minhash

import (
	"errors"
	"math"

	"github.com/will-rowe/ntHash"
)

// if set true, ntHash will return the canonical k-mer (inspects both strands for each k-mer and returns the lowest hash value)
const CANONICAL = false

/*
  The minHash struct contains all the minimum hash values for a sequence
*/
type minHash struct {
	kSize int
	signature []uint64
}

// Add a sequence to the minHash
func (minHash *minHash) Add(sequence []byte) error {
		// initiate the rolling hash
		hasher, err := ntHash.New(&sequence, minHash.kSize)
		if err != nil {
			return err
		}
	// get hashed kmers from read
		for hash := range hasher.Hash(CANONICAL) {
			// for each hashed k-mer, try adding it to the sketch
			for i, minVal := range minHash.signature {
				// split the hashed k-mer (uint64) into two uint32
				h1, h2 := uint32(hash), uint32(hash>>32)
				// get the new hash value for this signature position
				newVal := uint64(h1 + uint32(i)*h2)
				// evaluate and add to the signature if it is a minimum
				if newVal < minVal {
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
		return 0, errors.New("length of minhash signatures do not match\n")
	}
	intersect := 0
	for i := range minHash.signature {
		if minHash.signature[i] == querySig[i] {
			intersect++
		}
	}
	return float64(intersect) / float64(len(minHash.signature)), nil
}

// NewminHash initiates a minHash struct and populates the signature with max values
func NewMinHash(kSize, sigSize int) *minHash {
	signature := make([]uint64, sigSize)
	for i := range signature {
		signature[i] = math.MaxUint64
	}
	return &minHash{
		kSize: kSize,
		signature: signature,
	}
}
