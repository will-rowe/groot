package lshforest

import (
	"testing"
)

// test constructor
func TestConstructor(t *testing.T) {
	lshf := NewLSHforest(42, 0.99)
	funcs, buckets := lshf.Settings()
	if buckets != 1 || funcs != 42 {
		t.Fatal("incorrect number of buckets/funcs")
	}
}
