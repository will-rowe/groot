// testing is incomplete, more to be added...
package lshIndex

import (
    "fmt"
    "os"
	"testing"
)

var (
    // test graph windows
    entry1 = &GraphWindow{
		Key       : fmt.Sprintf("g%dn%do%d", 1, 2, 3),
		Size      : 100,
		Signature : []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
    }
    entry2 = &GraphWindow{
		Key       : fmt.Sprintf("g%dn%do%d", 1, 3, 1),
		Size      : 100,
		Signature : []uint64{1, 4, 3, 4, 5, 5, 7, 4, 9, 10},
    }
    entry3 = &GraphWindow{
		Key       : fmt.Sprintf("g%dn%do%d", 3, 22, 2),
		Size      : 100,
		Signature : []uint64{4, 4, 3, 4, 5, 6, 7, 4, 9, 4},
    }
    entries = []*GraphWindow{entry1, entry2, entry3}
    // LSH Forest parameters
    jsThresh = 0.85
    // LSH Ensemble parameters
    numPart = 4
    numHash = 10
    maxK = 4
    // query for LSH Forest
    query1 = &GraphWindow{
		Key       : fmt.Sprintf("g%dn%do%d", 1, 2, 3),
		Size      : 100,
		Signature : []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
    }
    // query for LSH Ensemble
    query2 = &GraphWindow{
		Key       : fmt.Sprintf("g%dn%do%d", 1, 2, 3),
		Size      : 50,
		Signature : []uint64{1, 1, 3, 4, 5, 6, 7, 8, 9, 10},
    }
)

// test the lshForest constructor, add a record and query it
func Test_lshForestConstructor(t *testing.T) {
    index := NewLSHforest(len(entry1.Signature), jsThresh)
    numHF, numBucks := index.Lshes[0].Settings()
	t.Logf("\tnumber of hash functions per bucket: %d\n", numHF)
    t.Logf("\tnumber of buckets: %d\n", numBucks)
    index.Add(entry1.Key, entry1.Signature, 0)
    index.Index()
    done := make(chan struct{})
    defer close(done)
    var check string
    for result := range index.Query(query1.Signature, query1.Size, jsThresh, done) {
        check = result.(string)
        if check != entry1.Key {
            t.Fatal()
        }
    }
    if check == "" {
        t.Fatal("no result from LSH Forest")
    }
}

// test the lshForest constructor and add a set of records, then query
func Test_lshForestBootstrap(t *testing.T) {
    index := NewLSHforest(len(entry1.Signature), jsThresh)
    for _, i := range entries {
        index.Add(i.Key, i.Signature, 0)
    }
    if len(index.Partitions) != 1 || index.SingleForest != true {
        t.Fatal()
    }
    index.Index()
    done := make(chan struct{})
    defer close(done)
    var check string
    for result := range index.Query(query1.Signature, query1.Size, jsThresh, done) {
        check = result.(string)
        if check != entry1.Key {
            t.Fatal("incorrect result returned from LSH Forest")
        }
    }
    if check == "" {
        t.Fatal("no result from LSH Forest")
    }
}

// test the lshForest dump and load methods
func Test_lshForestDump(t *testing.T) {
    index := NewLSHforest(len(entry1.Signature), jsThresh)
    for _, i := range entries {
        index.Add(i.Key, i.Signature, 0)
    }
    if err := index.Dump("./lsh.index"); err != nil {
        t.Fatal(err)
    }
    index2 := NewLSHforest(len(entry1.Signature), jsThresh)
    if err := index2.Load("./lsh.index"); err != nil {
        t.Fatal(err)
    }
    if err := os.Remove("./lsh.index"); err != nil {
        t.Fatal(err)
    }
    done := make(chan struct{})
    defer close(done)
    var check string
    for result := range index2.Query(query1.Signature, query1.Size, jsThresh, done) {
        check = result.(string)
        if check != entry1.Key {
            t.Fatal(check)
        }
    }
    if check == "" {
        t.Fatal("no result from LSH Forest")
    }
}


// test the lshEnsemble constructor, add the records and query it
func Test_lshEnsembleBootstrap(t *testing.T) {
    index := BootstrapLshEnsemble(numPart, numHash, maxK, len(entries), Windows2Chan(entries))
    index.Index()
    done := make(chan struct{})
    defer close(done)
    var check string
    for result := range index.Query(query2.Signature, query2.Size, jsThresh, done) {
        check = result.(string)
        if check != entry1.Key {
            t.Fatal("incorrect result returned from LSH Ensemble")
        }
    }
    if check == "" {
        t.Fatal("no result from LSH ensemble")
    }
}

// test the lshEnsemble dump and load methods
func Test_lshEnsembleDump(t *testing.T) {
    index := BootstrapLshEnsemble(numPart, numHash, maxK, len(entries), Windows2Chan(entries))
    if err := index.Dump("./lsh.index"); err != nil {
        t.Fatal(err)
    }
    index2 := NewLSHensemble(make([]Partition, numPart), numHash, maxK)
    if err := index2.Load("./lsh.index"); err != nil {
        t.Fatal(err)
    }
    if err := os.Remove("./lsh.index"); err != nil {
        t.Fatal(err)
    }
    done := make(chan struct{})
    defer close(done)
    var check string
    for result := range index2.Query(query2.Signature, query2.Size, jsThresh, done) {
        check = result.(string)
        if check != entry1.Key {
            t.Fatal()
        }
    }
    if check == "" {
        t.Fatal("no result from LSH ensemble")
    }
}
