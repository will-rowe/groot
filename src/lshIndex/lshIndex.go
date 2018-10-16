package lshIndex

import (
	//"errors"
	"github.com/orcaman/concurrent-map"
)

// set to 2/4/8 for 16bit/32bit/64bit hash values
const HASH_SIZE = 8
// integration precision for optimising number of bands + hash functions in LSH Forest
const PRECISION = 0.01
// number of partitions and maxK to use in LSH Ensemble (TODO: add these as customisable parameters for GROOT)
const PARTITIONS = 6
const MAXK = 4

// error messages
//var (
	//querySizeError = errors.New("Query size is > +/- 10 bases of reference windows, re-index using --containment")
//)

// NewLSHensemble initializes a new index consisting of MinHash LSH implemented using LshForest.
// numHash is the number of hash functions in MinHash.
// maxK is the maximum value for the MinHash parameter K - the number of hash functions per "band".
func NewLSHensemble(parts []Partition, numHash, maxK int) *LshEnsemble {
	lshes := make([]*LshForest, len(parts))
	for i := range lshes {
		lshes[i] = newLshForest(maxK, numHash/maxK)
	}
	return &LshEnsemble{
		Lshes:      lshes,
		Partitions: parts,
		MaxK:       maxK,
		NumHash:    numHash,
		paramCache: cmap.New(),
	}
}

// NewLshForest initializes a new index consisting of MinHash LSH implemented using a single LshForest.
// sigSize is the number of hash functions in MinHash.
// jsThresh is the minimum Jaccard similarity needed for a query to return a match
func NewLSHforest(sigSize int, jsThresh float64) *LshEnsemble {
	// calculate the optimal number of bands and hash functions to use
	numHashFuncs, numBands, _, _ := optimise(sigSize, jsThresh)
	lshes := make([]*LshForest, 1)
	lshes[0] = newLshForest(numHashFuncs, numBands)
	return &LshEnsemble{
		Lshes:      lshes,
		Partitions: make([]Partition, 1),
		MaxK:       numBands,
		NumHash:    numHashFuncs,
		paramCache: cmap.New(),
		SingleForest:	true,
	}
}

// BoostrapLshEnsemble builds an index from a channel of domains.
// The returned index consists of MinHash LSH implemented using LshForest.
// numPart is the number of partitions to create.
// numHash is the number of hash functions in MinHash.
// maxK is the maximum value for the MinHash parameter K - the number of hash functions per "band".
// GraphWindow is a channel emitting windows (don't need to be sorted by their sizes as windows are constant) TODO: should probably add a check for this
func BootstrapLshEnsemble(numPart, numHash, maxK, totalNumWindows int, windows <-chan *GraphWindow) *LshEnsemble {
	index := NewLSHensemble(make([]Partition, numPart), numHash, maxK)
	bootstrap(index, totalNumWindows, windows)
	return index
}

// bootstrap
func bootstrap(index *LshEnsemble, totalNumWindows int, windows <-chan *GraphWindow) {
	numPart := len(index.Partitions)
	depth := totalNumWindows / numPart
	var currDepth, currPart int
	for rec := range windows {
		index.Add(rec.Key, rec.Signature, currPart)
		currDepth++
		index.Partitions[currPart].Upper = rec.Size
		if currDepth >= depth && currPart < numPart-1 {
			currPart++
			index.Partitions[currPart].Lower = rec.Size
			currDepth = 0
		}
	}
	return
}

// Windows2Chan is a utility function that converts a GraphWindow slice in memory to a GraphWindow channel.
func Windows2Chan(windows []*GraphWindow) <-chan *GraphWindow {
	c := make(chan *GraphWindow, 1000)
	go func() {
		for _, w := range windows {
			c <- w
		}
		close(c)
	}()
	return c
}