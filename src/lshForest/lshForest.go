// Package lshforest is the indexing scheme used for GROOT
package lshforest

import (
	"encoding/binary"
	"fmt"
	"sort"
	"sync"
)

// HASH_SIZE is set to 2/4/8 for 16bit/32bit/64bit hash values
const HASH_SIZE = 8

// IndexWrapper is the base type for the GROOT index
// it provides access to an LSH forest index
type IndexWrapper struct {
	sketchSize     int                   // the size of the sketch that the LSH forest is indexing
	initHashTables []map[string][]string // the initial hash tables keep track of all SubSequences until the .index() method is called to populate the buckets and make the LSH forest searchable
	forest         *LSHforest            // the underlying LSH Forest
	sync.RWMutex                         // used to manage RW access to the keyLookup map
}

// NewLSHforest is the constructor function
func NewLSHforest(sketchSize int, jsThresh float64) *IndexWrapper {

	// calculate the optimal number of buckets and hash functions to use, based on the length of MinHash sketch and a Jaccard Similarity theshhold
	k, l, _, _ := optimise(sketchSize, jsThresh)

	// create the initial hash tables
	iht := make([]map[string][]string, l)
	for i := range iht {
		iht[i] = make(map[string][]string)
	}

	// create the protobuf compliant LSH forest
	forest := &LSHforest{
		K:         int32(k),
		L:         int32(l),
		KeyLookup: make(map[string]*Key),
	}

	// return the address of the LSH forest wrapper
	return &IndexWrapper{
		sketchSize:     sketchSize,
		initHashTables: iht,
		forest:         forest,
	}
}

// Settings will return the number of hash functions per partition and the total number of partitions used by the LSH Forest
func (IndexWrapper *IndexWrapper) Settings() (int32, int32) {
	return IndexWrapper.forest.K, IndexWrapper.forest.L
}

// Add a graph sketch to the LSH Forest
func (IndexWrapper *IndexWrapper) Add(graphWindow *Key) error {
	if len(graphWindow.Sketch) != IndexWrapper.sketchSize {
		return fmt.Errorf("cannot add sketch: wrong size for index")
	}

	// add the key to the lookup map, using the stringified version as the lookup
	// there may be multiple copies of the same window
	// - one graph+node+offset can have several subpaths to window
	// - or windows can be derived from identical regions of the graph that multiple sequences share
	stringifiedKey := fmt.Sprintf("g%dn%do%dp%d", graphWindow.GraphID, graphWindow.Node, graphWindow.OffSet, graphWindow.Ref)
	if err := IndexWrapper.addKey(stringifiedKey, graphWindow); err != nil {
		return err
	}

	// split the sketch into the right number of buckets, compressing each one down to a single string value
	partitionedSketch := make([]string, IndexWrapper.forest.L)
	for i := int32(0); i < IndexWrapper.forest.L; i++ {
		partitionedSketch[i] = CompressSketch2String(graphWindow.Sketch[i*IndexWrapper.forest.K : (i+1)*IndexWrapper.forest.K])
	}

	// iterate over each bucket of the initial hash tables in the LSH forest
	for i := 0; i < len(IndexWrapper.initHashTables); i++ {

		// if the current partition of the sketch isn't in the current bucket in the LSH forest, add it
		if _, ok := IndexWrapper.initHashTables[i][partitionedSketch[i]]; !ok {
			IndexWrapper.initHashTables[i][partitionedSketch[i]] = make([]string, 1)
			IndexWrapper.initHashTables[i][partitionedSketch[i]][0] = stringifiedKey
		} else {

			// if it is, append the current key (graph location) to this hashed sketch bucket
			IndexWrapper.initHashTables[i][partitionedSketch[i]] = append(IndexWrapper.initHashTables[i][partitionedSketch[i]], stringifiedKey)
		}
	}
	return nil
}

// GetKey will return the Key for the stringified version
func (IndexWrapper *IndexWrapper) GetKey(keystring string) (*Key, error) {
	IndexWrapper.Lock()
	returnKey, ok := IndexWrapper.forest.KeyLookup[keystring]
	IndexWrapper.Unlock()
	if ok {
		return returnKey, nil
	}
	return nil, fmt.Errorf("key not found in LSH Forest: %v", keystring)
}

// addKey will add a Key to the lookup map, storing it under a stringified version of the key
func (IndexWrapper *IndexWrapper) addKey(stringifiedKey string, key *Key) error {
	if stringifiedKey == "" {
		return fmt.Errorf("empty key passed to the LSH Forest")
	}
	IndexWrapper.Lock()
	if _, contains := IndexWrapper.forest.KeyLookup[stringifiedKey]; contains {
		return fmt.Errorf("duplicated key encountered: %v", stringifiedKey)
	}
	IndexWrapper.forest.KeyLookup[stringifiedKey] = key
	IndexWrapper.Unlock()
	return nil
}

// clearSketches deletes all the sketches stored with each key in the KeyLookup - this is to save some space after they are no longer needed
func (IndexWrapper *IndexWrapper) clearSketches() {
	for _, key := range IndexWrapper.forest.KeyLookup {
		key.Sketch = []uint64{}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// bucket is a partition in the LSHforest
// it contains the subsequences (and the keys) for an index region (i.e. region i..K for sketch A, B, C)
// Pair is defined in lshforest.proto
type bucket []*Pair

// methods to satisfy the sort interface, allowing each bucket to be searchable
func (bucket bucket) Len() int           { return len(bucket) }
func (bucket bucket) Swap(i, j int)      { bucket[i], bucket[j] = bucket[j], bucket[i] }
func (bucket bucket) Less(i, j int) bool { return bucket[i].SubSequence < bucket[j].SubSequence }

/*
// index transfers the contents of each initialHashTable to the bucket arrays so they can be sorted and searched
func (IndexWrapper *IndexWrapper) index() {

	// make sure the LSH Forest is empty
	IndexWrapper.forest.Buckets = make([]*Bucket, IndexWrapper.forest.L)
	for i := range IndexWrapper.forest.Buckets {
		IndexWrapper.forest.Buckets[i] = &Bucket{}
	}

	// gather the initial hash tables, put them into tmp buckets, sort them and then add them to the LSH Forest
	for i := 0; i < len(IndexWrapper.initHashTables); i++ {

		// empty the hash table into the tmp bucket
		tmpBucket := make(bucket, len(IndexWrapper.initHashTables[i]))
		counter := 0
		for seq, keys := range IndexWrapper.initHashTables[i] {
			tmpBucket[counter] = &Pair{SubSequence: seq, Keys: keys}
			counter++
		}
		sort.Sort(tmpBucket)

		// convert it to the protobuff compliant Bucket
		bucket := &Bucket{}
		bucket.Pairs = tmpBucket

		// add the sorted bucket to the LSH Forest
		IndexWrapper.forest.Buckets[i] = bucket

		// clear the initial hashtable that has just been processed
		IndexWrapper.initHashTables[i] = nil
	}
}

*/

// transferToBuckets transfers the contents of each initialHashTable to the bucket arrays so they can be written to disk, sorted and searched
func (IndexWrapper *IndexWrapper) transferToBuckets() {

	// make sure the LSH Forest is empty
	IndexWrapper.forest.Buckets = make([]*Bucket, IndexWrapper.forest.L)
	for i := range IndexWrapper.forest.Buckets {
		IndexWrapper.forest.Buckets[i] = &Bucket{}
	}

	// gather the initial hash tables, put them into tmp buckets, sort them and then add them to the LSH Forest
	for i := int32(0); i < int32(len(IndexWrapper.initHashTables)); i++ {

		// empty the hash table into the tmp bucket
		tmpBucket := make(bucket, len(IndexWrapper.initHashTables[i]))
		counter := 0
		for compressedBucket, stringifiedKeys := range IndexWrapper.initHashTables[i] {
			key, err := IndexWrapper.GetKey(stringifiedKeys[0])
			if err != nil {
				panic(err)
			}
			tmpBucket[counter] = &Pair{SubSequence: compressedBucket, Keys: stringifiedKeys, SketchPartition: key.Sketch[i*IndexWrapper.forest.K : (i+1)*IndexWrapper.forest.K]}
			counter++
		}

		// sort the bucket based on the compressedBucket (string representation of the sketch partition)
		sort.Sort(tmpBucket)

		// convert it to the protobuff compliant Bucket
		for _, tmp := range tmpBucket {

			// protobuff can't store string representations of binary arrays - results in UTF-8 error
			// we have to re-calculate this on index load.... TODO: this is garbage, is there another way?
			tmp.SubSequence = ""
		}
		bucket := &Bucket{}
		bucket.Pairs = tmpBucket

		// add the sorted bucket to the LSH Forest
		IndexWrapper.forest.Buckets[i] = bucket

		// clear the initial hashtable that has just been processed
		IndexWrapper.initHashTables[i] = nil
	}

	// clear the sketches in the lookup map as they are no longer needed
	IndexWrapper.clearSketches()
}

// convertSketchPartitionsToString takes each []uint64 in each bucket and compresses it from []uint64 -> binary -> string
func (IndexWrapper *IndexWrapper) convertSketchPartitionsToString() {

	// TODO: check that transferToBuckets() has been run

	// for each bucket, convert each entry's SketchPartition to a stringified representation
	for _, bucket := range IndexWrapper.forest.Buckets {
		for i := 0; i < len(bucket.Pairs); i++ {
			bucket.Pairs[i].SubSequence = CompressSketch2String(bucket.Pairs[i].SketchPartition)

			// remove the sketch partition as we only need the stringified version now
			bucket.Pairs[i].SketchPartition = []uint64{}
		}
	}
}

// CompressSketch2String is a function to compress an array of unsigned integers as a string - taken from https://github.com/ekzhu/minhash-lsh
func CompressSketch2String(sketch []uint64) string {
	s := make([]byte, HASH_SIZE*len(sketch))
	buf := make([]byte, 8)
	for i, v := range sketch {
		binary.LittleEndian.PutUint64(buf, v)
		copy(s[i*HASH_SIZE:(i+1)*HASH_SIZE], buf[:HASH_SIZE])
	}
	return string(s)
}
