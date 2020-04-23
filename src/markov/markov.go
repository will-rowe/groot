package markov

import (
	"errors"
	"fmt"
	"sync"
)

// Tokens are wrapped around a sequence of words to maintain the start and end transition counts
const (
	StartToken = "$"
	EndToken   = "^"
)

// Chain is a markov chain instance
type Chain struct {
	Order        int
	statePool    *spool
	frequencyMat map[int]sparseArray
	lock         *sync.RWMutex
}

// NewChain creates an instance of Chain
func NewChain(order int) *Chain {
	chain := Chain{Order: order}
	chain.statePool = &spool{
		stringMap: make(map[string]int),
		intMap:    make(map[int]string),
	}
	chain.frequencyMat = make(map[int]sparseArray, 0)
	chain.lock = new(sync.RWMutex)
	return &chain
}

// Add adds the transition counts to the chain for a given sequence of words
func (chain *Chain) Add(input []string) {
	startTokens := array(StartToken, chain.Order)
	endTokens := array(EndToken, chain.Order)
	tokens := []string{}
	tokens = append(tokens, startTokens...)
	tokens = append(tokens, input...)
	tokens = append(tokens, endTokens...)

	pairs := MakePairs(tokens, chain.Order)
	for i := 0; i < len(pairs); i++ {
		pair := pairs[i]
		currentIndex := chain.statePool.add(pair.CurrentState.key())
		nextIndex := chain.statePool.add(pair.NextState)
		chain.lock.Lock()
		if chain.frequencyMat[currentIndex] == nil {
			chain.frequencyMat[currentIndex] = make(sparseArray, 0)
		}
		chain.frequencyMat[currentIndex][nextIndex]++

		chain.lock.Unlock()
	}
}

// Scale will adjust the transition count for a given N-gram + string transition count
func (chain *Chain) Scale(next string, current NGram, scaling float64) error {
	if scaling == 0.0 {
		return nil
	}
	startToken := true
	for i := 0; i < len(current); i++ {
		if current[i] == "" {
			if startToken {
				current[i] = StartToken
			} else {
				current[i] = EndToken
			}
		} else {
			startToken = false
		}
	}
	currentIndex := chain.statePool.add(current.key())
	nextIndex := chain.statePool.add(next)
	chain.lock.Lock()
	if chain.frequencyMat[currentIndex] == nil {
		return fmt.Errorf("current node has not yet been added to the markov chain")
	}
	chain.frequencyMat[currentIndex][nextIndex] *= scaling
	if chain.frequencyMat[currentIndex][nextIndex] < 0 {
		chain.frequencyMat[currentIndex][nextIndex] = 0
	}
	chain.lock.Unlock()
	return nil
}

// TransitionProbability returns the transition probability between two states
func (chain *Chain) TransitionProbability(next string, current NGram) (float64, error) {
	if len(current) != chain.Order {
		return 0, errors.New("N-gram length does not match chain order")
	}
	currentIndex, currentExists := chain.statePool.get(current.key())
	nextIndex, nextExists := chain.statePool.get(next)
	if !currentExists || !nextExists {
		return 0, nil
	}
	arr := chain.frequencyMat[currentIndex]
	sum := float64(arr.sum())
	freq := float64(arr[nextIndex])
	return freq / sum, nil
}
