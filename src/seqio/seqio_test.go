package seqio

import (
	"testing"
)

// setup variables
var (
	l1         = []byte("@0_chr1_0_186027_186126_263_(Bla)BIC-1:GQ260093:1-885:885")
	l2         = []byte("acagcaggaaggcttactggagaaacgtatcgactataagaatcgggtgatggaacctcactctcccatcagcgcacaacatagttcgacgggtatgacc")
	l3         = []byte("+")
	l4         = []byte("====@==@AAD?>D@@==DACBC?@BB@C==AB==A@D>AD==?CB==@=B?=A>D?=DB=?>>D@EB===??=@C=?C>@>@B>=?C@@>=====?@>=")
	kmerSize   = 7
	sketchSize = 10
)

// test results
var (
	expectedUpperCase  = []byte("ACAGCAGGAAGGCTTACTGGAGAAACGTATCGACTATAAGAATCGGGTGATGGAACCTCACTCTCCCATCAGCGCACAACATAGTTCGACGGGTATGACC")
	expectedTrimmedSeq = []byte("GAAGGCTTACTGGAGAAACGTATCGACTATAAGAATCGGGTGATGGAACCTCACTCTCCCATCAGCGCACAACATAGTTCGAC")
	expectedRevComp    = []byte("GTCGAACTATGTTGTGCGCTGATGGGAGAGTGAGGTTCCATCACCCGATTCTTATAGTCGATACGTTTCTCCAGTAAGCCTTC")
)

// test functions to check equality of slices
func ByteSliceCheck(a, b []byte) bool {
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
func Uint64SliceCheck(a, b []uint64) bool {
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// begin the tests
func TestReadConstructor(t *testing.T) {
	_, err := NewFASTQread(l1, l2, l3, l4)
	if err != nil {
		t.Fatalf("could not generate FASTQ read using NewFASTQread")
	}
}

func TestSeqMethods(t *testing.T) {
	read, err := NewFASTQread(l1, l2, l3, l4)
	if err != nil {
		t.Fatalf("could not generate FASTQ read using NewFASTQread")
	}
	read.BaseCheck()
	if ByteSliceCheck(read.Seq, expectedUpperCase) == false {
		t.Errorf("Bases2Upper method failed")
	}
	read.QualTrim(30)
	if ByteSliceCheck(read.Seq, expectedTrimmedSeq) == false {
		t.Errorf("QualTrim method failed")
	}
	read.RevComplement()
	if ByteSliceCheck(read.Seq, expectedRevComp) == false {
		t.Errorf("RevComplement method failed")
	}
}

func TestMinHash(t *testing.T) {
	read, err := NewFASTQread(l1, l2, l3, l4)
	if err != nil {
		t.Fatalf("could not generate FASTQ read using NewFASTQread")
	}
	// sketch using KHF MinHash
	if sketch, err := read.RunMinHash(kmerSize, sketchSize, true, nil); err != nil {
		t.Fatal(err)
	} else {
		t.Log(sketch)
	}
	// sketch using KMV MinHash
	if sketch, err := read.RunMinHash(kmerSize, sketchSize, false, nil); err != nil {
		t.Fatal(err)
	} else {
		t.Log(sketch)
	}
}
