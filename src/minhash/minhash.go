// Package minhash contains implementations of bottom-k and kmv MinHash algorithms. These implementations use the ntHash rolling hash function.
package minhash

// CANONICAL tell ntHash to return the canonical k-mer (this is used in the KMV sketch)
const CANONICAL bool = true

// MinHash is an interface to group the different flavours of MinHash implemented here
type MinHash interface {
	AddSequence([]byte) error
	GetSketch() []uint64
}

// GetReadSketch is a function to sketch a read sequence
func GetReadSketch(seq []byte, kmerSize, sketchSize uint, kmv bool) ([]uint64, error) {

	// create the MinHash data structure, using the specified algorithm flavour
	var mh MinHash
	if kmv {
		mh = NewKMVsketch(uint(kmerSize), uint(sketchSize))
	} else {
		mh = NewKHFsketch(uint(kmerSize), uint(sketchSize))
	}

	// use the AddSequence method to populate the MinHash
	err := mh.AddSequence(seq)

	// get the sketch
	sketch := mh.GetSketch()

	// if the sketch isn't at capacity (in the case of BottomK sketches), fill up the remainder with 0s
	if kmv && len(sketch) != int(sketchSize) {
		padding := make([]uint64, int(sketchSize)-len(sketch))
		for i := 0; i < len(padding); i++ {
			padding[i] = 0
		}
		sketch = append(sketch, padding...)

	}

	// return the MinHash sketch and any error
	return sketch, err
}

// seqNT4table is used to convert "ACGTN" to 01234 - from minimap2
var seqNT4table = [256]uint8{
	0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
}

// hash64 is a hash function for uint64 encoded k-mers (lifted from minimap2)
func hash64(key, mask uint64) uint64 {
	key = (^key + (key << 21)) & mask
	key = key ^ key>>24
	key = ((key + (key << 3)) + (key << 8)) & mask
	key = key ^ key>>14
	key = ((key + (key << 2)) + (key << 4)) & mask
	key = key ^ key>>28
	key = (key + (key << 31)) & mask
	return key
}

// splitmix64 is a 64-bit finalizer, used here as a second hash func for uint64 endcoded k-mers
func splitmix64(key uint64) uint64 {
	key = (key ^ (key >> 31) ^ (key >> 62)) * uint64(0x319642b2d24d8ec3)
	key = (key ^ (key >> 27) ^ (key >> 54)) * uint64(0x96de1b173f119089)
	key = key ^ (key >> 30) ^ (key >> 60)
	return key
}
