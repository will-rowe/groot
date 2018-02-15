#!/bin/bash

# install the software
go build

# index the example
./groot index -p 1 -i paper/data/clustered-ARG-databases/arg-annot-90 -o test-index

# align the test reads
./groot align -p 1 -i test-index -f testing/full-argannot-perfect-reads-small.fq.gz > out.bam

## TODO: check the output bam file is as expected
