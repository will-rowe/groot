// Package minhash contains implementations of bottom-k and kmv MinHash algorithms. These implementations use the nthash rolling hash function.
package minhash

// CANONICAL tell nthash to return the canonical k-mer (this is used in the KMV sketch)
const CANONICAL bool = true

// MinHash is an interface to group the different flavours of MinHash implemented here
type MinHash interface {
	AddSequence([]byte) error
	GetSketch() []uint64
}
