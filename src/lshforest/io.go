package lshforest

import (
	"fmt"
	"io/ioutil"

	proto "github.com/golang/protobuf/proto"
)

// Dump is a method to dump the LSHforest to file
func (IndexWrapper *IndexWrapper) Dump(path string) error {

	// converts the initial hash tables to sorted arrays
	IndexWrapper.transferToBuckets()

	// marshall the data
	data, err := proto.Marshal(IndexWrapper.forest)
	if err != nil {
		return err
	}

	// write to disk, returning any error
	return ioutil.WriteFile(path, data, 0644)
}

// Load is a method to load LSHforest from file
func (IndexWrapper *IndexWrapper) Load(path string) error {
	data, err := ioutil.ReadFile(path)
	if err != nil {
		return err
	}
	return IndexWrapper.LoadFromBytes(data)
}

// LoadFromBytes is a method to populate an LSH Forest instance using a byte slice
func (IndexWrapper *IndexWrapper) LoadFromBytes(data []byte) error {
	if len(data) == 0 {
		return fmt.Errorf("no data received to load the LSH Forest from")
	}

	// get the import holder ready
	importedForest := &LSHforest{}

	// unmarshall into the holder
	if err := proto.Unmarshal(data, importedForest); err != nil {
		return err
	}

	// some quick sanity checks
	if (importedForest.K * importedForest.L) != int32(IndexWrapper.sketchSize) {
		return fmt.Errorf("groot graphs and lsh forest were not created in same run")
	}
	if len(importedForest.Buckets) < 1 {
		return fmt.Errorf("LSH Forest is corrupted")
	}

	// attach the forest to the wrapper
	IndexWrapper.forest = importedForest

	// re-stringify all the sketch partitions in the LSH Forest buckets
	IndexWrapper.convertSketchPartitionsToString()
	return nil
}
