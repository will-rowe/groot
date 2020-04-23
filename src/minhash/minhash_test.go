package minhash

import (
	"testing"

	"github.com/adam-hanna/arrayOperations"
)

var (
	hashvalues      = []uint64{12345, 54321, 9999999, 98765}
	kmerSize        = uint(7)
	sketchSize      = uint(10)
	seqA            = []byte("ACTGCGTGCGTGAAACGTGCACGTGACGTG")
	seqArcomplement = []byte("CACGTCACGTGCACGTTTCACGCACGCAGT")
)

// kmerShredder is a helper function for yielding k-mers from a sequence
func kmerShredder(seq []byte, k uint) []string {
	numKmers := len(seq) - int(k) + 1
	kmers := make([]string, numKmers)
	for i := 0; i < numKmers; i++ {
		kmers[i] = string(seq[i : i+int(k)])
	}
	return kmers
}

// calJS is a helper function to calculate the jaccard similarity of two slices of strings
func calJS(setA, setB []string) float64 {
	if len(setA) != len(setB) {
		panic("set size mismatch")
	}
	intersect := intersection(setA, setB)
	return float64(len(intersect)) / float64(len(setA))
}

// intersection returns the common elements between two slices of strings
func intersection(a, b []string) []string {
	z, ok := arrayOperations.Intersect(a, b)
	if !ok {
		panic("Cannot find intersect")
	}

	slice, ok := z.Interface().([]string)
	if !ok {
		panic("Cannot convert to slice")
	}
	return slice
}

// BloomFilter test
func TestBloomfilter(t *testing.T) {
	filter := NewBloomFilter(3)
	for i := 0; i < len(hashvalues); i++ {
		filter.Add(hashvalues[i])
	}
	for i := 0; i < len(hashvalues); i++ {
		if !filter.Check(hashvalues[i]) {
			t.Fatalf("'%d' should be have been marked present", hashvalues[i])
		}
	}
	filter.Reset()
	for i := 0; i < len(hashvalues); i++ {
		if filter.Check(hashvalues[i]) {
			t.Fatalf("'%d' shouldn't be marked as present", hashvalues[i])
		}
	}
}

// Constructor test
func TestMinHashConstructors(t *testing.T) {
	mhKHF := NewKHFsketch(kmerSize, sketchSize)
	if len(mhKHF.GetSketch()) != int(sketchSize) || mhKHF.sketchSize != sketchSize || mhKHF.kmerSize != kmerSize {
		t.Fatalf("NewKHFsketch constructor did not initiate MinHash KHF sketch correctly")
	}
	mhKMV := NewKMVsketch(kmerSize, sketchSize)
	if mhKMV.sketchSize != sketchSize || mhKMV.kmerSize != kmerSize {
		t.Fatalf("NewKMVsketch constructor did not initiate MinHash KMV sketch correctly")
	}
}

// Add test for KHF
func TestKHFadd(t *testing.T) {
	mhKHF := NewKHFsketch(kmerSize, sketchSize)

	// try adding a sequence that is too short for the given k
	if err := mhKHF.AddSequence(seqA[0:1]); err == nil {
		t.Fatal("should fault as sequences must be >= kmerSize")
	}

	// try adding a sequence that passes the length check
	if err := mhKHF.AddSequence(seqA); err != nil {
		t.Fatal(err)
	}
}

// Add test for KMV
func TestKMVadd(t *testing.T) {
	mhKMV := NewKMVsketch(kmerSize, sketchSize)

	// try adding a sequence that is too short for the given k
	if err := mhKMV.AddSequence(seqA[0:1]); err == nil {
		t.Fatal("should fault as sequences must be >= kmerSize")
	}

	// try adding a sequence that passes the length check
	if err := mhKMV.AddSequence(seqA); err != nil {
		t.Fatal(err)
	}
}

func TestSimilarityEstimates(t *testing.T) {

	// get actual JS for the test sequences
	setA := kmerShredder(seqA, kmerSize)
	setB := kmerShredder(seqArcomplement, kmerSize)
	js := calJS(setA, setB)
	t.Logf("actual Jaccard similarity: %.2f", js)

	// test KHF
	mhKHF1 := NewKHFsketch(kmerSize, sketchSize)
	if err := mhKHF1.AddSequence(seqA); err != nil {
		t.Fatal(err)
	}
	mhKHF2 := NewKHFsketch(kmerSize, sketchSize)
	if err := mhKHF2.AddSequence(seqArcomplement); err != nil {
		t.Fatal(err)
	}

	// test KMV
	mhKMV1 := NewKMVsketch(kmerSize, sketchSize)
	if err := mhKMV1.AddSequence(seqA); err != nil {
		t.Fatal(err)
	}
	mhKMV2 := NewKMVsketch(kmerSize, sketchSize)
	if err := mhKMV2.AddSequence(seqArcomplement); err != nil {
		t.Fatal(err)
	}

	// because these sketches used canonical k-mers, and the two sequences being compared are reverse complements of each other, they should yield identical k-mer sets
	// for KMV, the sketch size == the set size, which means Jaccard is guarenteed to be 1.0, if this has been implemented right! This guarentee doesn't hold for KHF though
	khfJS, err := mhKHF1.GetSimilarity(mhKHF2)
	if err != nil {
		t.Fatal(err)
	}
	if khfJS != 1.0 {
		t.Fatalf("similarity estimate from KHF MinHash should be 1.0, not: %.2f", khfJS)
	}

	kmvJS, err := mhKMV1.GetSimilarity(mhKMV2)
	if err != nil {
		t.Fatal(err)
	}
	if kmvJS != 1.0 {
		t.Fatalf("similarity estimate from KMV MinHash should be 1.0, not: %.2f", kmvJS)
	}

}

// benchmark KHF
func BenchmarkKHF(b *testing.B) {
	mhKHF1 := NewKHFsketch(kmerSize, sketchSize)

	// run the add method b.N times
	for n := 0; n < b.N; n++ {
		if err := mhKHF1.AddSequence(seqA); err != nil {
			b.Fatal(err)
		}
	}
}

// benchmark KMV
func BenchmarkKMV(b *testing.B) {
	mhKMV1 := NewKMVsketch(kmerSize, sketchSize)

	// run the add method b.N times
	for n := 0; n < b.N; n++ {
		if err := mhKMV1.AddSequence(seqA); err != nil {
			b.Fatal(err)
		}
	}
}
