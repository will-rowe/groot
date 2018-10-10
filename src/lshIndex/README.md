# lshIndex package

Since version 0.8.0, GROOT has two options for indexing the variation graphs - lshForest or lshEnsemble. With the addition of lshEnsemble, GROOT can now receive variable read lengths and seed these against graphs using containment search.

This has required a significant re-write of the lshForest code, which is all contained in this directory. The majority of this code now comes from the [lshEnsemble package](https://godoc.org/github.com/ekzhu/lshensemble) by ekzhu. I have just made a few changes:

* removed methods that were unnecessary for GROOT
* added a method to write the index to disk
* added methods to generate a single LSH Forest index using a Jaccard Similarity threshold and signature length parameter
