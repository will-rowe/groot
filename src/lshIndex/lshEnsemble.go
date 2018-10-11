package lshIndex

import (
	"encoding/gob"
	"fmt"
	"os"
	"sync"

	"github.com/orcaman/concurrent-map"
	"github.com/will-rowe/groot/src/seqio"
)

type param struct {
	k int
	l int
}

// Partition represents a domain size partition in the LSH Ensemble index.
type Partition struct {
	Lower int `json:"lower"`
	Upper int `json:"upper"`
}

// KeyLookupMap relates the stringified seqio.Key to the original, allowing LSH index search results to easily be related to graph locations
type KeyLookupMap map[string]seqio.Key

// GraphWindow represents a region of a variation graph
type GraphWindow struct {
	// The unique key of this window
	Key interface{}
	// The window size
	Size int
	// The MinHash signature of this window
	Signature []uint64
}

// LshEnsemble represents an LSH Ensemble index.
type LshEnsemble struct {
	Partitions []Partition
	Lshes      []*LshForest
	MaxK       int
	NumHash    int
	paramCache cmap.ConcurrentMap
	Indexed bool
	SingleForest bool
	KeyLookup   KeyLookupMap
}

// Add a new domain to the index given its partition ID - the index of the partition.
// The added domain won't be searchable until the Index() function is called.
func (e *LshEnsemble) Add(key interface{}, sig []uint64, partInd int) {
	e.Lshes[partInd].Add(key, sig)
}

// Index makes all added domains searchable.
func (e *LshEnsemble) Index() {
	for i := range e.Lshes {
		e.Lshes[i].Index()
	}
	e.Indexed = true
}

// Query returns the candidate domain keys in a channel.
// This function is given the MinHash signature of the query domain, sig, the domain size,
// the containment threshold, and a cancellation channel.
// Closing channel done will cancel the query execution.
// The query signature must be generated using the same seed as the signatures of the indexed domains,
// and have the same number of hash functions.
func (e *LshEnsemble) Query(sig []uint64, size int, threshold float64, done <-chan struct{}) <-chan interface{} {
	if e.SingleForest {
		return e.queryForest(sig, done)
	} 
	params := e.computeParams(size, threshold)
	return e.queryWithParam(sig, params, done)
}

//
func (e *LshEnsemble) queryWithParam(sig []uint64, params []param, done <-chan struct{}) <-chan interface{} {
	// Collect candidates from all partitions
	keyChan := make(chan interface{})
	var wg sync.WaitGroup
	wg.Add(len(e.Lshes))
	for i := range e.Lshes {
		go func(lsh *LshForest, k, l int) {
			lsh.Query(sig, k, l, keyChan, done)
			wg.Done()
		}(e.Lshes[i], params[i].k, params[i].l)
	}
	go func() {
		wg.Wait()
		close(keyChan)
	}()
	return keyChan
}

//
func (e *LshEnsemble) queryForest(sig []uint64, done <-chan struct{}) <-chan interface{} {
	keyChan := make(chan interface{})
	var wg sync.WaitGroup
	wg.Add(1)
	go func(lsh *LshForest) {
		lsh.Query(sig, -1, -1, keyChan, done)
		wg.Done()
	}(e.Lshes[0])
	go func() {
		wg.Wait()
		close(keyChan)
	}()
	return keyChan
}

// Compute the optimal k and l for each partition
func (e *LshEnsemble) computeParams(size int, threshold float64) []param {
	params := make([]param, len(e.Partitions))
	for i, p := range e.Partitions {
		x := p.Upper
		key := cacheKey(x, size, threshold)
		if cached, exist := e.paramCache.Get(key); exist {
			params[i] = cached.(param)
		} else {
			optK, optL, _, _ := e.Lshes[i].OptimalKL(x, size, threshold)
			computed := param{optK, optL}
			e.paramCache.Set(key, computed)
			params[i] = computed
		}
	}
	return params
}

// Make a cache key with threshold precision to 2 decimal points
func cacheKey(x, q int, t float64) string {
	return fmt.Sprintf("%.8x %.8x %.2f", x, q, t)
}

// Dump an LSH index to disk
func (LshEnsemble *LshEnsemble) Dump(path string) error {
	if LshEnsemble.Indexed == true {
		return fmt.Errorf("cannot dump the LSH Index after running the indexing method")
	}
	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()
	encoder := gob.NewEncoder(file)
	if err := encoder.Encode(LshEnsemble); err != nil {
		return err
	 }
	for _, lsh := range LshEnsemble.Lshes {
		for _, bandContents := range lsh.initHashTables {
			err := encoder.Encode(bandContents)
			if err != nil {
				return err
			}
		}
	}
	err = encoder.Encode(LshEnsemble.KeyLookup)
	if err != nil {
		return err
	}
	return nil
}

// Load an LSH index from disk
func (LshEnsemble *LshEnsemble) Load(path string) error {
	file, err := os.Open(path)
	if err != nil {
		return err
	}
	defer file.Close()
	decoder := gob.NewDecoder(file)
	if err := decoder.Decode(&LshEnsemble); err != nil {
		return(err)
	 }
	for _, lsh := range LshEnsemble.Lshes {
		for _, bandContents := range lsh.initHashTables {
			err := decoder.Decode(&bandContents)
			if err != nil {
				return err
			}
		}
	}
	err = decoder.Decode(&LshEnsemble.KeyLookup)
	if err != nil {
		return err
	}
	LshEnsemble.Index()
	return nil
}
