// Package graph is used to process graphs. It converts, writes, aligns reads and processes GROOT graphs.
package graph

import (
	"bytes"
	"encoding/gob"
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"sort"
	"sync"

	"github.com/ekzhu/lshensemble"
	"github.com/will-rowe/groot/src/lshforest"
)

// ContainmentIndex is a wrapper for the LSH Ensemble data structure
type ContainmentIndex struct {

	// LookupMap relates the stringified keys in the index to the graph windows
	LookupMap map[string]*lshforest.Key

	// DomainRecords is used for construction of the index, then destroyed
	DomainRecords []*lshensemble.DomainRecord

	// NumPart is the number of partitions
	NumPart int

	// MaxK is the number of hash funcs per band
	MaxK int

	// WindowSize is the number of k-mers in the graph windows
	WindowSize int

	// SketchSize is the size of the sketches being indexed (num hash funcs)
	SketchSize int

	// LSHensemble is the LSH Ensemble index
	LSHensemble *lshensemble.LshEnsemble

	// unexported:
	numSketches int          // number of sketches in the containment index
	lock        sync.RWMutex // lock is a mutex to manage RW access to the LookupMap
}

/*

TODO:

add the Key struct here and remove LSH Forest guff


*/

// PrepareIndex prepares a LSH Ensemble index from a map of domain records and a map of keys-windows
func PrepareIndex(domainRecMap map[int]*lshensemble.DomainRecord, lookupMap map[string]*lshforest.Key, numPart, maxK, windowSize, sketchSize int) *ContainmentIndex {

	// setup the struct
	ci := &ContainmentIndex{
		DomainRecords: make([]*lshensemble.DomainRecord, len(domainRecMap)),
		LookupMap:     lookupMap,
		NumPart:       numPart,
		MaxK:          maxK,
		WindowSize:    windowSize,
		SketchSize:    sketchSize,
	}

	// get the domain records into a sorted arrays
	for i, rec := range domainRecMap {
		ci.DomainRecords[i] = rec
	}
	sort.Sort(lshensemble.BySize(ci.DomainRecords))
	return ci
}

// Dump is a method to write a containment index to disk
func (ContainmentIndex *ContainmentIndex) Dump(filePath string) error {

	// make sure it has had the prepare method run
	if ContainmentIndex.DomainRecords == nil || ContainmentIndex.LookupMap == nil {
		return fmt.Errorf("must run PrepareIndex before dumping index to disk")
	}

	// error if the index has already been loaded once
	if ContainmentIndex.numSketches != 0 {
		return fmt.Errorf("this index cannot be dumped a second time")
	}

	// write to disk
	fh, err := os.Create(filePath)
	defer fh.Close()
	if err != nil {
		return err
	}
	encoder := gob.NewEncoder(fh)
	return encoder.Encode(ContainmentIndex)
}

// Load is a method to load a containment index from disk and populate the LSH Ensemble
func (ContainmentIndex *ContainmentIndex) Load(filePath string) error {
	data, err := ioutil.ReadFile(filePath)
	if err != nil {
		return err
	}
	if len(data) == 0 {
		return fmt.Errorf("index appears empty")
	}

	return ContainmentIndex.LoadFromBytes(data)
}

// LoadFromBytes is a method to load the containment index from a byte array
func (ContainmentIndex *ContainmentIndex) LoadFromBytes(data []byte) error {
	buf := bytes.NewBuffer(data)
	decoder := gob.NewDecoder(buf)
	var err error
	if err := decoder.Decode(ContainmentIndex); err != nil {
		return err
	}
	if ContainmentIndex.DomainRecords == nil {
		return fmt.Errorf("loaded an empty index file")
	}

	/*
		// populate the LSH Ensemble
		ContainmentIndex.LSHensemble, err = lshensemble.BootstrapLshEnsembleOptimal(ContainmentIndex.NumPart, ContainmentIndex.SketchSize, ContainmentIndex.MaxK,
			func() <-chan *lshensemble.DomainRecord {
				return lshensemble.Recs2Chan(ContainmentIndex.DomainRecords)
			})
	*/

	// Create index using equi-depth partitioning
	// You can also use BootstrapLshEnsemblePlusEquiDepth for better accuracy
	ContainmentIndex.LSHensemble, err = lshensemble.BootstrapLshEnsembleEquiDepth(ContainmentIndex.NumPart, ContainmentIndex.SketchSize, ContainmentIndex.MaxK, len(ContainmentIndex.DomainRecords), recs2Chan(ContainmentIndex.DomainRecords))

	ContainmentIndex.numSketches = len(ContainmentIndex.DomainRecords)

	// get rid of the domain records as they are not needed anymore
	ContainmentIndex.DomainRecords = nil

	return err
}

// Query is temp function to check the the index can be queried
func (ContainmentIndex *ContainmentIndex) Query(querySig []uint64, querySize int, containmentThreshold float64) ([]*lshforest.Key, error) {
	done := make(chan struct{})
	defer close(done)
	results := []*lshforest.Key{}
	for hit := range ContainmentIndex.LSHensemble.Query(querySig, querySize, containmentThreshold, done) {

		key, err := ContainmentIndex.getKey(hit.(string))
		if err != nil {
			return nil, err
		}

		// full containment check
		// TODO: this should be optional
		if lshensemble.Containment(querySig, key.Sketch, querySize, ContainmentIndex.WindowSize) > containmentThreshold {
			results = append(results, key)
		}
	}
	return results, nil
}

// getKey will return the Key for the stringified version
func (ContainmentIndex *ContainmentIndex) getKey(keystring string) (*lshforest.Key, error) {
	ContainmentIndex.lock.Lock()
	returnKey, ok := ContainmentIndex.LookupMap[keystring]
	ContainmentIndex.lock.Unlock()
	if ok {
		return returnKey, nil
	}
	return nil, fmt.Errorf("key not found in LSH Ensemble: %v", keystring)
}

// recs2Chan is a utility function that converts a DomainRecord slice in memory to a DomainRecord channel
func recs2Chan(recs []*lshensemble.DomainRecord) <-chan *lshensemble.DomainRecord {
	c := make(chan *lshensemble.DomainRecord, 1000)
	go func() {
		for i, r := range recs {

			// this is a massive kludge to fix the out of memory issue on index load in wasm
			// see: https://github.com/golang/go/issues/32840#issuecomment-506929883
			if runtime.GOOS == "js" && i%10000 == 0 {
				runtime.GC()
			}
			c <- r
		}
		close(c)
	}()
	return c
}
