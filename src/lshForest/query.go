package lshforest

import (
	"sort"
	"sync"
)

// Query is the exported method for querying and returning similar sketches from the LSH forest
func (IndexWrapper *IndexWrapper) Query(sketch []uint64) []string {
	result := make([]string, 0)

	// more info on done chans for explicit cancellation in concurrent pipelines: https://blog.golang.org/pipelines
	done := make(chan struct{})

	// collect query results and aggregate in a single array to send back
	for key := range IndexWrapper.runQuery(sketch, done) {
		result = append(result, key)
	}
	close(done)
	return result
}

// runQuery does the actual work
func (IndexWrapper *IndexWrapper) runQuery(sketch []uint64, done <-chan struct{}) <-chan string {
	queryResultChan := make(chan string)
	go func() {
		defer close(queryResultChan)

		//  convert the query sketch from []uint64 to a string
		stringifiedSketch := make([]string, IndexWrapper.forest.L)
		for i := int32(0); i < IndexWrapper.forest.L; i++ {
			stringifiedSketch[i] = CompressSketch2String(sketch[i*IndexWrapper.forest.K : (i+1)*IndexWrapper.forest.K])
		}

		// don't send back multiple copies of the same key
		seens := make(map[string]bool)

		// compress internal nodes using a prefix
		prefixSize := HASH_SIZE * (IndexWrapper.forest.K - 1)

		// run concurrent hashtable queries
		keyChan := make(chan string)
		var wg sync.WaitGroup
		wg.Add(int(IndexWrapper.forest.L))
		for i := int32(0); i < IndexWrapper.forest.L; i++ {
			go func(bucket []*Pair, queryChunk string) {
				defer wg.Done()

				// sort.Search uses binary search to find and return the smallest index i in [0, n) at which f(i) is true
				index := sort.Search(len(bucket), func(x int) bool { return bucket[x].SubSequence[:prefixSize] >= queryChunk[:prefixSize] })

				// k is the index returned by the search
				if index < len(bucket) && bucket[index].SubSequence[:prefixSize] == queryChunk[:prefixSize] {
					for j := index; j < len(bucket) && bucket[j].SubSequence[:prefixSize] == queryChunk[:prefixSize]; j++ {
						if bucket[j].SubSequence == queryChunk {

							// if the query matches the bucket, send the keys as search results
							for _, key := range bucket[j].Keys {
								select {
								case keyChan <- key:
								case <-done:
									return
								}
							}
						}
					}
				}
			}(IndexWrapper.forest.Buckets[i].Pairs, stringifiedSketch[i])
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
