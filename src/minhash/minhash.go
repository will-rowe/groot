// the minhash package contains a MinHash implementation that is adapted from go-minhash (https://godoc.org/github.com/dgryski/go-minhash)
package minhash

import (
	"errors"
	"math"
)

/*
  The MinHash struct contains all the minimum hash values for a sequence
*/
type MinHash struct {
	signature []uint64
	hashFunc1 Hash64
	hashFunc2 Hash64
}
type Hash64 func([]byte) uint64

/*
  A method to hash a k-mer and add signature to the MinHash struct
*/
func (self *MinHash) Add(kmer []byte) {
	hashValue1, hashValue2 := self.hashFunc1(kmer), self.hashFunc2(kmer)
	for i, minVal := range self.signature {
		newVal := hashValue1 + uint64(i)*hashValue2
		if newVal < minVal {
			self.signature[i] = newVal
		}
	}
}

/*
  A method to dump the MinHash signature
*/
func (self *MinHash) Signature() []uint64 {
	return self.signature
}

/*
  A method to estimate Jaccard Similarity between a MinHash struct and a query MH signature
*/
func (self *MinHash) Similarity(querySig []uint64) (float64, error) {
	if len(self.signature) != len(querySig) {
		return 0, errors.New("length of minhash signatures do not match\n")
	}
	intersect := 0
	for i := range self.signature {
		if self.signature[i] == querySig[i] {
			intersect++
		}
	}
	return float64(intersect) / float64(len(self.signature)), nil
}

/*
  A function to create a new MinHash struct
*/
func NewMinHash(h1, h2 Hash64, size int) *MinHash {
	signature := make([]uint64, size)
	for i := range signature {
		signature[i] = math.MaxUint64
	}
	newMinHash := new(MinHash)
	newMinHash.hashFunc1 = h1
	newMinHash.hashFunc2 = h2
	newMinHash.signature = signature
	return newMinHash
}
