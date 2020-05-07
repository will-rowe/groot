// Package lshe is used to index the graphs. It wraps the LSH ensemble library: https://github.com/ekzhu/lshensemble
package lshe

import (
	"bytes"
	"encoding/gob"
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"sync"

	"github.com/ekzhu/lshensemble"
)

// Key relates sketches of reads and graph traversals to specific windows of a graph
type Key struct {
	GraphID        uint32             // identifies the graph
	Node           uint64             // identifies the first node in the graph window
	OffSet         uint32             // identifies the offset of a window within the first node
	ContainedNodes map[uint64]float64 // describes the traversal through the graph for the window
	Ref            []uint32           // the IDs for the reference sequences that contains this window
	RC             bool               // identifies if the read has been reverse complemented (NOT USED)
	Sketch         []uint64           // the sketch of this graph window
	Freq           float64            // records the number of k-mers this graph window has received during read mapping
	MergeSpan      uint32             // indicates maximum distance between graph windows this key represents (used in window merging if sketches identical)
	WindowSize     uint32             // the size of the window that was sketched (prior to merging)
}

// Keys is used to hold multiple keys (and satisfies the sort interface)
type Keys []Key

func (k Keys) Len() int           { return len(k) }
func (k Keys) Less(i, j int) bool { return k[i].Node < k[j].Node }
func (k Keys) Swap(i, j int)      { k[i], k[j] = k[j], k[i] }

// ContainmentIndex is a wrapper for the LSH Ensemble data structure
type ContainmentIndex struct {
	NumPart        int            // NumPart is the number of partitions
	MaxK           int            // MaxK is the number of hash funcs per band
	NumWindowKmers int            // NumWindowKmers is the number of k-mers in the graph windows
	SketchSize     int            // SketchSize is the size of the sketches being indexed (num hash funcs)
	WindowLookup   map[string]Key // WindowLookup is a map linking windows to a their sketch in the LSH Ensemble

	// unexported
	lshEnsemble *lshensemble.LshEnsemble // lshEnsemble is the LSH Ensemble index
	numSketches int                      // numSketches in the containment index
	locker      sync.Mutex               // locker to control access to the WindowLookup
}

// InitIndex will get a containment index struct ready
func InitIndex(numPart, maxK, numWindowKmers, sketchSize int) *ContainmentIndex {
	return &ContainmentIndex{
		NumPart:        numPart,
		MaxK:           maxK,
		NumWindowKmers: numWindowKmers,
		SketchSize:     sketchSize,
		WindowLookup:   make(map[string]Key),
	}
}

// AddWindow will add a window to a ContainmentIndex
func (ContainmentIndex *ContainmentIndex) AddWindow(windowLookup string, window Key) error {
	if _, ok := ContainmentIndex.WindowLookup[windowLookup]; ok {
		return fmt.Errorf("duplicate window key can't be inserted into index: %v", windowLookup)
	}
	ContainmentIndex.WindowLookup[windowLookup] = window
	return nil
}

// Dump is a method to write a containment index to disk
func (ContainmentIndex *ContainmentIndex) Dump(filePath string) error {

	// make sure it has had the prepare method run
	if len(ContainmentIndex.WindowLookup) == 0 {
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
	ContainmentIndex.numSketches = len(ContainmentIndex.WindowLookup)
	if ContainmentIndex.numSketches == 0 {
		return fmt.Errorf("loaded an empty index file")
	}

	// convert the sketched windows into domain records and a lookup map - put recs into a channel for LSH ensemble bootstrapping
	// NOTE: not sorting the domain records as all were built using same window size
	recChan := make(chan *lshensemble.DomainRecord, 100)
	go func() {
		i := 0
		for windowStringKey, window := range ContainmentIndex.WindowLookup {

			// this is a massive kludge to fix the out of memory issue on index load in wasm
			// see: https://github.com/golang/go/issues/32840#issuecomment-506929883
			if runtime.GOOS == "js" && i%10000 == 0 {
				runtime.GC()
			}

			// send the domain record for this window
			recChan <- &lshensemble.DomainRecord{
				windowStringKey,
				ContainmentIndex.NumWindowKmers,
				window.Sketch,
			}
			i++
		}
		close(recChan)
	}()

	// create index using equi-depth partitioning
	ContainmentIndex.lshEnsemble, err = lshensemble.BootstrapLshEnsembleEquiDepth(ContainmentIndex.NumPart, ContainmentIndex.SketchSize, ContainmentIndex.MaxK, ContainmentIndex.numSketches, recChan)

	return err
}

// Query wraps the LSH ensemble query method
// query sig is the sketch
// query size is the number of k-mers in the query sequence
// containment threshold is the containment threshold...
func (ContainmentIndex *ContainmentIndex) Query(querySig []uint64, querySize int, containmentThreshold float64) (map[uint32]Keys, error) {
	done := make(chan struct{})
	defer close(done)
	results := make(map[uint32]Keys)
	for hit := range ContainmentIndex.lshEnsemble.Query(querySig, querySize, containmentThreshold, done) {
		key, err := ContainmentIndex.getKey(hit.(string))
		if err != nil {
			return nil, err
		}

		// full containment check
		// TODO: this should probably be optional but overhead seems minimal
		if lshensemble.Containment(querySig, key.Sketch, querySize, ContainmentIndex.NumWindowKmers) > containmentThreshold {
			if len(results[key.GraphID]) == 0 {
				results[key.GraphID] = Keys{key}
			} else {
				results[key.GraphID] = append(results[key.GraphID], key)
			}
		}
	}

	return results, nil
}

// getKey will return the Key for the stringified version
func (ContainmentIndex *ContainmentIndex) getKey(keystring string) (Key, error) {
	ContainmentIndex.locker.Lock()
	returnKey, ok := ContainmentIndex.WindowLookup[keystring]
	ContainmentIndex.locker.Unlock()
	if ok {
		return returnKey, nil
	}
	return Key{}, fmt.Errorf("key not found in LSH Ensemble: %v", keystring)
}
