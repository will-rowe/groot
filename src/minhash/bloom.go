package minhash

import "sync"

// defaultSize used to create a bloom filter
const defaultSize = 10000

var mask [64]uint64

// init will prepare the mask prior to creating a bloom filter
func init() {
	mask[0] = 1
	for i := 1; i < len(mask); i++ {
		mask[i] = 2 * mask[i-1]
	}
}

// BloomFilter is the bloom filter type
type BloomFilter struct {
	size   uint64
	sketch []uint64
	lock   sync.RWMutex // lock the node for read/write access
}

// Reset will clear all marked bits in the Bloom Filter sketch
func (BloomFilter *BloomFilter) Reset() {
	for i := 0; i < len(BloomFilter.sketch); i++ {
		BloomFilter.sketch[i] = 0
	}
}

// Add is a method to add a hashed k-mer to the Bloom Filter sketch
func (BloomFilter *BloomFilter) Add(kmer uint64) {
	h := (kmer % BloomFilter.size)
	c := h / 64 // cell
	o := h % 64 // offset
	BloomFilter.lock.Lock()
	BloomFilter.sketch[c] = BloomFilter.sketch[c] | mask[o]
	BloomFilter.lock.Unlock()
}

// Check is a method to check a hashed k-mer against the Bloom Filter sketch
func (BloomFilter *BloomFilter) Check(kmer uint64) bool {
	h := (kmer % BloomFilter.size)
	c := h / 64 // cell
	o := h % 64 // offset
	BloomFilter.lock.Lock()
	defer BloomFilter.lock.Unlock()
	return (BloomFilter.sketch[c] & mask[o]) > 0
}

// NewBloomFilter is a Bloom Filter constructor, using a specified sized
func NewBloomFilter(size int) *BloomFilter {
	if size > 64 {
		size = size / 64
	} else {
		size = 1
	}
	return &BloomFilter{
		size:   64 * uint64(size),
		sketch: make([]uint64, size),
	}
}

// NewDefaultBloomFilter is a Bloom Filter constructor, using the default size
func NewDefaultBloomFilter() *BloomFilter {
	return NewBloomFilter(defaultSize)
}
