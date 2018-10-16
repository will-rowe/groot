package minhash

import (
	"testing"
)

var (
	kSize = 11
	sigSize = 24
	sequence = []byte("ACTGCGTGCGTGAAACGTGCACGTGACGTG")
)

func TestMinHashConstructor(t *testing.T) {
	mh := NewMinHash(kSize, sigSize)
	if len(mh.Signature()) != sigSize {
		t.Fatalf("NewMinHash did not initiate correctly")
	}
}

func TestAdd(t *testing.T) {
	mh := NewMinHash(kSize, sigSize)
	// try adding a sequence that is too short for the given k
	if err := mh.Add(sequence[0:3]); err == nil {
		t.Fatal("should fault as sequences must be >= kSize")
	}
	// try adding a sequence that passes the legnth check
	err := mh.Add(sequence)
	if err != nil {
		t.Fatal(err)
	}
}

func TestSimilarity(t *testing.T) {
	mh := NewMinHash(kSize, sigSize)
	_ = mh.Add(sequence)
	mh2 := NewMinHash(kSize, sigSize)
	_ = mh2.Add(sequence)
	// make sure JS calculation works
	js, err := mh.Similarity(mh2.Signature())
	if err != nil {
		t.Fatal(err)
	}
	if js != 1.0 {
		t.Fatal("incorrect JS calculation")
	}
	// make sure the method checks work
	mh3 := NewMinHash(kSize, (sigSize+1))
	_ = mh.Add(sequence)
	_, err = mh.Similarity(mh3.Signature())
	if err == nil {
		t.Fatal("should fault as signature lengths vary")
	}
}