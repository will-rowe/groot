package minhash

import (
	"container/heap"
	"fmt"
	"sort"

	"github.com/will-rowe/ntHash"
)

// KMVsketch is the structure for the K-Minimum Values MinHash sketch of a set of k-mers
type KMVsketch struct {
	kmerSize   uint
	sketchSize uint
	sketch     []uint64
	heap       *IntHeap
}

// NewKMVsketch is the constructor for a KMVsketch data structure
func NewKMVsketch(k, s uint) *KMVsketch {
	newSketch := &KMVsketch{
		kmerSize:   k,
		sketchSize: s,
		heap:       &IntHeap{},
	}

	// init the heap
	heap.Init(newSketch.heap)
	return newSketch
}

// AddSequence is a method to decompose a read to canonical kmers, hash them and add any minimums to the sketch
func (KMVsketch *KMVsketch) AddSequence(sequence []byte) error {

	// check the sequence length
	if len(sequence) < int(KMVsketch.kmerSize) {
		return fmt.Errorf("sequence length (%d) is short than k-mer length (%d)", len(sequence), KMVsketch.kmerSize)
	}

	// initiate the rolling ntHash
	hasher, err := ntHash.New(&sequence, KMVsketch.kmerSize)
	if err != nil {
		return err
	}

	// get hashed kmers from sequence and evaluate
	for hv := range hasher.Hash(CANONICAL) {

		// if the heap isn't full yet, go ahead and add the hash
		if len(*KMVsketch.heap) < int(KMVsketch.sketchSize) {
			heap.Push(KMVsketch.heap, hv)

			// re-establish the heap ordering after adding the new hash
			heap.Fix(KMVsketch.heap, 0)

			// or if the incoming hash is smaller than the hash at the top of the heap, add the hash and remove the larger one from the heap
		} else if hv < (*KMVsketch.heap)[0] {

			// replace the largest value currently in the sketch with the new hash
			(*KMVsketch.heap)[0] = hv

			// re-establish the heap ordering after adding the new hash
			heap.Fix(KMVsketch.heap, 0)
		}
	}
	return nil
}

// GetSketch is a method to convert the heap to a []uint64 sketch, which has been sorted (max > min)
func (KMVsketch *KMVsketch) GetSketch() []uint64 {
	sketch := make(IntHeap, len(*KMVsketch.heap))
	copy(sketch, *KMVsketch.heap)
	sort.Sort(sketch)
	return sketch
}

// GetSimilarity estimates the similarity between two k-mer sets based on the KMV sketch
func (mh1 *KMVsketch) GetSimilarity(mh2 MinHash) (float64, error) {

	// check this is a pair of KMV
	if fmt.Sprintf("%T", mh1) != fmt.Sprintf("%T", mh2) {
		return 0.0, fmt.Errorf("mismatched MinHash types: %T vs. %T", mh1, mh2)
	}
	querySketch, ok := mh2.(*KMVsketch)
	if !ok {
		return 0.0, fmt.Errorf("could not assert sketch is a KMV")
	}

	if len(*mh1.heap) != len(*querySketch.heap) {
		panic("sketches do not have the same number of minimums")
	}

	mins := make(map[uint64]int, len(*mh1.heap))
	for _, v := range *mh1.heap {
		mins[v]++
	}
	intersect := 0

	for _, v := range *querySketch.heap {
		if count, ok := mins[v]; ok && count > 0 {
			intersect++
			mins[v] = count - 1
		}
	}

	maxlength := len(*mh1.heap)
	if maxlength < len(*querySketch.heap) {
		maxlength = len(*querySketch.heap)
	}

	return (float64(intersect) / float64(maxlength)), nil
}
