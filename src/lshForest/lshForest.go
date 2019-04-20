package lshForest

import (
	"encoding/binary"
	"fmt"
	"io/ioutil"
	"math"
	"sort"
	"sync"

	"github.com/will-rowe/groot/src/seqio"
	"gopkg.in/vmihailenco/msgpack.v2"
)

// set to 2/4/8 for 16bit/32bit/64bit hash values
const HASH_SIZE = 8

/*
	The LSH forest
*/
type LSHforest struct {
	K                   int
	L                   int
	InitialHashTables   []initialHashTable
	hashTables          []hashTables
	hashedSignatureFunc hashedSignatureFunc
	KeyLookup           KeyLookupMap
}

/*
	A function to construct the LSH forest
*/
func NewLSHforest(sigSize int, jsThresh float64) *LSHforest {
	// calculate the optimal number of bands and hash functions to use based on the length of MinHash signature and a Jaccard Similarity theshhold
	K, L, _, _ := optimise(sigSize, jsThresh)
	// create the initial hash tables
	InitialHashTables := make([]initialHashTable, L)
	for i := range InitialHashTables {
		InitialHashTables[i] = make(initialHashTable)
	}
	// create the hash tables that will be populated once the LSH forest indexing method has been run
	indexedHashTables := make([]hashTables, L)
	for i := range indexedHashTables {
		indexedHashTables[i] = make(hashTables, 0)
	}
	// create the KeyLookup map to relate signatures to graph locations
	KeyLookup := make(KeyLookupMap)
	// return the address of the new LSH forest
	newLSHforest := new(LSHforest)
	newLSHforest.K = K
	newLSHforest.L = L
	newLSHforest.InitialHashTables = InitialHashTables
	newLSHforest.hashTables = indexedHashTables
	newLSHforest.hashedSignatureFunc = hashedSignatureFuncGen(HASH_SIZE)
	newLSHforest.KeyLookup = KeyLookup
	return newLSHforest
}

/*
	The types needed by the LSH forest
*/
// this map relates the stringified seqio.Key to the original, allowing LSHforest search results to easily be related to graph locations
type KeyLookupMap map[string]seqio.Key

// graphKeys is a slice containing all the stringified graphKeys for a given hashed signature
type graphKeys []string

// the initial hash table uses the hashed signature as a key - the values are the corresponding graphKeys
type initialHashTable map[string]graphKeys

// a band is a single hash table that is stored in the indexedHashTables - it contains the band of a hash signature and the corresponding graphKeys
type band struct {
	HashedSignature string
	graphKeys       graphKeys
}

// this is populated during indexing -- it is a slice of bands and can be sorted
type hashTables []band

//methods to satisfy the sort interface
func (h hashTables) Len() int           { return len(h) }
func (h hashTables) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h hashTables) Less(i, j int) bool { return h[i].HashedSignature < h[j].HashedSignature }

// the hashkey function type and the generator function
type hashedSignatureFunc func([]uint64) string

func hashedSignatureFuncGen(hashValueSize int) hashedSignatureFunc {
	return func(sig []uint64) string {
		hashedSig := make([]byte, hashValueSize*len(sig))
		buf := make([]byte, 8)
		for i, v := range sig {
			// use the ByteOrder interface to write binary data
			// use the LittleEndian implementation and call the Put method
			binary.LittleEndian.PutUint64(buf, v)
			copy(hashedSig[i*hashValueSize:(i+1)*hashValueSize], buf[:hashValueSize])
		}
		return string(hashedSig)
	}
}

/*
	A method to return the number of hash functions and number of bands set by the LSH forest
*/
func (self *LSHforest) Settings() (K, L int) {
	return self.K, self.L
}

/*
	A method to add a minhash signature and graph key to the LSH forest
*/
func (self *LSHforest) Add(key string, sig []uint64) {
	// split the signature into the right number of bands and then hash each one
	hashedSignature := make([]string, self.L)
	for i := 0; i < self.L; i++ {
		hashedSignature[i] = self.hashedSignatureFunc(sig[i*self.K : (i+1)*self.K])
	}
	// iterate over each band in the LSH forest
	for i := 0; i < len(self.InitialHashTables); i++ {
		// if the current band in the signature isn't in the current band in the LSH forest, add it
		if _, ok := self.InitialHashTables[i][hashedSignature[i]]; !ok {
			self.InitialHashTables[i][hashedSignature[i]] = make(graphKeys, 1)
			self.InitialHashTables[i][hashedSignature[i]][0] = key
			// if it is, append the current key (graph location) to this hashed signature band
		} else {
			self.InitialHashTables[i][hashedSignature[i]] = append(self.InitialHashTables[i][hashedSignature[i]], key)
		}
	}
}

/*
	A method to index the graph (transfers contents of each initialHashTable so they can be sorted and searched)
*/
func (self *LSHforest) Index() {
	// iterate over the empty indexed hash tables
	for i := range self.hashTables {
		// transfer contents from the corresponding band in the initial hash table
		for HashedSignature, keys := range self.InitialHashTables[i] {
			self.hashTables[i] = append(self.hashTables[i], band{HashedSignature, keys})
		}
		// sort the new hashtable and store it in the corresponding slot in the indexed hash tables
		sort.Sort(self.hashTables[i])
		// clear the initial hashtable that has just been processed
		self.InitialHashTables[i] = make(initialHashTable)
	}
}

/*
	Methods to dump the LSH forest to disk and then load it again
*/
// Dump an LSH index to disk
func (self *LSHforest) Dump(path string) error {
	if len(self.hashTables[0]) != 0 {
		return fmt.Errorf("cannot dump the LSH Forest after running the indexing method")
	}
	b, err := msgpack.Marshal(self)
	if err != nil {
		return err
	}
	return ioutil.WriteFile(path, b, 0644)
}

// Load an LSH index from disk
func (self *LSHforest) Load(path string) error {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		return err
	}
	return msgpack.Unmarshal(b, self)
}

/*
	A method to query a MinHash signature against the LSH forest
*/
func (self *LSHforest) Query(sig []uint64) []string {
	result := make([]string, 0)
	// more info on done chans for explicit cancellation in concurrent pipelines: https://blog.golang.org/pipelines
	done := make(chan struct{})
	defer close(done)
	// collect query results and aggregate in a single array to send back
	for key := range self.runQuery(sig, done) {
		result = append(result, key)
	}
	return result
}

func (self *LSHforest) runQuery(sig []uint64, done <-chan struct{}) <-chan string {
	queryResultChan := make(chan string)
	go func() {
		defer close(queryResultChan)
		//  hash the query signature
		hashedSignature := make([]string, self.L)
		for i := 0; i < self.L; i++ {
			hashedSignature[i] = self.hashedSignatureFunc(sig[i*self.K : (i+1)*self.K])
		}
		// don't send back multiple copies of the same key
		seens := make(map[string]bool)
		// compress internal nodes using a prefix
		prefixSize := HASH_SIZE * self.K
		// run concurrent hashtable queries
		keyChan := make(chan string)
		var wg sync.WaitGroup
		wg.Add(self.L)
		for i := 0; i < self.L; i++ {
			go func(band hashTables, queryChunk string) {
				defer wg.Done()
				// sort.Search uses binary search to find and return the smallest index i in [0, n) at which f(i) is true
				index := sort.Search(len(band), func(x int) bool { return band[x].HashedSignature[:prefixSize] >= queryChunk })
				// k is the index returned by the search
				if index < len(band) && band[index].HashedSignature[:prefixSize] == queryChunk {
					for j := index; j < len(band) && band[j].HashedSignature[:prefixSize] == queryChunk; j++ {
						// copies key values from this hashtable to the keyChan until all values from band[j] copied or done is closed
						for _, key := range band[j].graphKeys {
							select {
							case keyChan <- key:
							case <-done:
								return
							}
						}
					}
				}
			}(self.hashTables[i], hashedSignature[i])
		}
		go func() {
			wg.Wait()
			close(keyChan)
		}()
		for key := range keyChan {
			if _, seen := seens[key]; seen {
				continue
			}
			queryResultChan <- key
			seens[key] = true
		}
	}()
	return queryResultChan
}

//  the following funcs are taken from https://github.com/ekzhu/minhash-lsh

// optimise returns the optimal number of hash functions and the optimal number of bands for Jaccard similarity search, as well as  the false positive and negative probabilities.
func optimise(sigSize int, jsThresh float64) (int, int, float64, float64) {
	optimumK, optimumL := 0, 0
	fp, fn := 0.0, 0.0
	minError := math.MaxFloat64
	for l := 1; l <= sigSize; l++ {
		for k := 1; k <= sigSize; k++ {
			if l*k > sigSize {
				break
			}
			currFp := probFalsePositive(l, k, jsThresh, 0.01)
			currFn := probFalseNegative(l, k, jsThresh, 0.01)
			currErr := currFn + currFp
			if minError > currErr {
				minError = currErr
				optimumK = k
				optimumL = l
				fp = currFp
				fn = currFn
			}
		}
	}
	return optimumK, optimumL, fp, fn
}

// Compute the integral of function f, lower limit a, upper limit l, and
// precision defined as the quantize step
func integral(f func(float64) float64, a, b, precision float64) float64 {
	var area float64
	for x := a; x < b; x += precision {
		area += f(x+0.5*precision) * precision
	}
	return area
}

// Probability density function for false positive
func falsePositive(l, k int) func(float64) float64 {
	return func(j float64) float64 {
		return 1.0 - math.Pow(1.0-math.Pow(j, float64(k)), float64(l))
	}
}

// Probability density function for false negative
func falseNegative(l, k int) func(float64) float64 {
	return func(j float64) float64 {
		return 1.0 - (1.0 - math.Pow(1.0-math.Pow(j, float64(k)), float64(l)))
	}
}

// Compute the cummulative probability of false negative given threshold t
func probFalseNegative(l, k int, t, precision float64) float64 {
	return integral(falseNegative(l, k), t, 1.0, precision)
}

// Compute the cummulative probability of false positive given threshold t
func probFalsePositive(l, k int, t, precision float64) float64 {
	return integral(falsePositive(l, k), 0, t, precision)
}
