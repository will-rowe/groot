package minhash

// IntHeap is a min-heap of uint64s (we're satisfying the heap interface: https://golang.org/pkg/container/heap/)
type IntHeap []uint64

// the less method is returning the larger value, so that it is at index position 0 in the heap
func (IntHeap IntHeap) Less(i, j int) bool { return IntHeap[i] > IntHeap[j] }
func (IntHeap IntHeap) Swap(i, j int)      { IntHeap[i], IntHeap[j] = IntHeap[j], IntHeap[i] }
func (IntHeap IntHeap) Len() int           { return len(IntHeap) }

// Push is a method to add an element to the heap
func (IntHeap *IntHeap) Push(x interface{}) {
	// dereference the pointer to modify the slice's length, not just its contents
	*IntHeap = append(*IntHeap, x.(uint64))
}

// Pop is a method to remove an element from the heap
func (IntHeap *IntHeap) Pop() interface{} {
	old := *IntHeap
	n := len(old)
	x := old[n-1]
	*IntHeap = old[0 : n-1]
	return x
}
