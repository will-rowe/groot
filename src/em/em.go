// Package em is the groot implementation of the expectation-maximization algorithm for finding the most likely paths through the graphs
// Currently, there is no weighting for each path, as they are all so similar in length and the nodes have already been weighted. This may change after I've thought on it some more...
package em

import (
	"fmt"
	"math"
)

// tolerance
var tolerance = math.Nextafter(1, 2) - 1

// EMrunner is the type to run the EM algorithm
type EMrunner struct {
	paths   map[uint32][]byte   // the paths in the graph (i.e. the possible reference genes)
	lengths map[uint32]int      // the length of each path
	ecMap   map[uint64][]uint32 // equivalence class (ec) map -> key=ecID, values=nodeIDs
	counts  map[uint64]float64  // the count for each ec, where index == ec ID
	//weights []float64 // the weight for each ec
	numIterations int // the max number of iterations for the EM
	minIterations int // the minimum number of iterations for the EM
	alpha         []float64
	alphaB4zeroes []float64
	rho           []float64
	iterationsRan int // counts how many iterations were run once Run() is called
}

// NewEM is the constructor
func NewEM(nIterations int, mIterations int, paths map[uint32][]byte, lengths map[uint32]int, ecMap map[uint64][]uint32, counts map[uint64]float64) (*EMrunner, error) {
	if nIterations < mIterations {
		return nil, fmt.Errorf("number of EM iterations (%d) must be greater than minimum iterations (%d)", nIterations, mIterations)
	}
	numPaths := len(paths)
	// set up alpha
	ud := 1.0 / float64(numPaths)
	alpha := make([]float64, numPaths)
	alphaB4zeroes := make([]float64, numPaths)
	rho := make([]float64, numPaths)
	for i := 0; i < numPaths; i++ {
		alpha[i] = ud
		alphaB4zeroes[i] = ud
		rho[i] = 0.0
	}

	// set up the EM
	return &EMrunner{
		paths:         paths,
		ecMap:         ecMap,
		counts:        counts,
		lengths:       lengths,
		numIterations: nIterations,
		minIterations: mIterations,
		alpha:         alpha,
		alphaB4zeroes: alphaB4zeroes,
		rho:           rho,
	}, nil
}

// Run is the method to run the EM algorithm
func (EMrunner *EMrunner) Run() error {
	numPaths := len(EMrunner.paths)
	nextAlpha := make([]float64, numPaths)
	for i := 0; i < numPaths; i++ {
		nextAlpha[i] = 0.0
	}
	denom := 0.0
	alphaLimit := 1e-7
	alphaChange := 1e-2
	alphaChangeLimit := 1e-2
	finalRound := false

	// start the EM
	numEMiterations := 0
	for numEMiterations = 0; numEMiterations < EMrunner.numIterations; numEMiterations++ {

		// iterate over the equivalence classes
		for ec, v := range EMrunner.ecMap {
			denom = 0.0
			if _, ok := EMrunner.counts[ec]; !ok {
				return fmt.Errorf("could no look up count for EC")
			}
			if EMrunner.counts[ec] == 0 {
				continue
			}

			// get the number of genes in this ec
			numEC := len(v)

			// get the weight
			//wv := EMrunner.weights[ec]

			// compute the dominator
			for pathIterator := 0; pathIterator < numEC; pathIterator++ {
				//denom += EMrunner.alpha[v[pathIterator]] * wv[pathIterator]
				denom += EMrunner.alpha[v[pathIterator]]
			}
			if denom < tolerance {
				continue
			}

			// compute the update step
			countNorm := EMrunner.counts[ec] / denom
			for pathIterator := 0; pathIterator < numEC; pathIterator++ {
				//nextAlpha[v[pathIterator]] += (wv[pathIterator] * EMrunner.alpha[v[pathIterator]]) * countNorm
				nextAlpha[v[pathIterator]] += EMrunner.alpha[v[pathIterator]] * countNorm
			}
		}

		stopEM := false
		chcount := 0
		for ec := 0; ec < numPaths; ec++ {
			if nextAlpha[ec] > alphaChangeLimit && (math.Abs(nextAlpha[ec]-EMrunner.alpha[ec])/nextAlpha[ec]) > alphaChange {
				chcount++
			}

			// reassign alpha_ to next_alpha
			EMrunner.alpha[ec] = nextAlpha[ec]

			// clear all next_alpha values 0 for next iteration
			nextAlpha[ec] = 0.0
		}
		if chcount == 0 && numEMiterations > EMrunner.minIterations {
			stopEM = true
		}
		if finalRound {
			break
		}
		if stopEM {
			finalRound = true
			EMrunner.alphaB4zeroes = EMrunner.alphaB4zeroes[:len(EMrunner.alpha)]
			for ec := 0; ec < numPaths; ec++ {
				EMrunner.alphaB4zeroes[ec] = EMrunner.alpha[ec]
				if EMrunner.alpha[ec] < alphaLimit/10.0 {
					EMrunner.alpha[ec] = 0.0
				}
			}
		}
	} // end of EM

	// if EM ran for the maximum number of iterations
	if numEMiterations == EMrunner.numIterations {
		EMrunner.alphaB4zeroes = EMrunner.alphaB4zeroes[:len(EMrunner.alpha)]
		for ec := 0; ec < numPaths; ec++ {
			EMrunner.alphaB4zeroes[ec] = EMrunner.alpha[ec]
		}
	}
	EMrunner.iterationsRan = numEMiterations

	return nil
}

// Return is the method to return EM results to the haplotype pipeline process
func (EMrunner *EMrunner) Return() (int, []float64, error) {
	if EMrunner.iterationsRan < 1 {
		return 0, nil, fmt.Errorf("no EM iterations were ran")
	}
	return EMrunner.iterationsRan, EMrunner.alpha, nil
}
