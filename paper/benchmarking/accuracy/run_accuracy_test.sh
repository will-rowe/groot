#!/bin/env bash
#
# This test checks the seeding and alignment accuracy of the program
#
# The CAMI metagenomic dataset is spiked with ARGs from the CARD database
# Randomreads are generated at the specified coverage level (errors are allowed)
# Groot and ARGOAP are both run
#
# run with parallel using a command like:
# cat test-read-set-sizes | parallel --gnu "sh run_accuracy_test.sh {}"
#
# combine the usage outputs into a single tsv:
# printf "CPU usage:\ttime (wall clock seconds)\n" > usage-data.tsv
# cat usage-for* >> usage-data.tsv && rm usage-for*
#
# Requires:
#  BBMap
#  Go 1.9
#  gtime (brew install gnu-time)
#

# test parameters
THREADS=8
READ_LEN=150
K_SIZE=7
SIG_SIZE=128
JS=0.99
NUM_READS=$1
echo "test parameters:"
echo "  - k-mer size: " $K_SIZE
echo "  - signature size: " $SIG_SIZE
echo "  - Jaccard similarity: " $JS

mkdir ${NUM_READS} && cd ${NUM_READS}

# create the test reads
echo "creating test reads..."
randomreads.sh ref=../../../data/full-ARG-databases/arg-annot-db/argannot-args.fna out=${NUM_READS}-test-reads.fq length=$READ_LEN  reads=$NUM_READS maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 adderrors=false > /dev/null 2>&1

# install GROOT
echo "installing GROOT..."
go build ../../../..
go build ../accuracy-test.go

# index the ARGannot database
echo "indexing the ARG-annot database..."
gtime -f "\tmax. resident set size (kb): %M\n\tCPU usage: %P\n\ttime (wall clock): %E\n\ttime (CPU seconds): %S\n" ./groot index -i ../../../data/clustered-ARG-databases/arg-annot-90 -o groot-index -l $READ_LEN -k $K_SIZE -s $SIG_SIZE -j $JS -p $THREADS

# align the test reads
echo "aligning reads..."
gtime -o  ../usage-for-${NUM_READS}-reads.tsv -f "%P\t%e\n" ./groot align -i groot-index -f ${NUM_READS}-test-reads.fq -p $THREADS > groot.bam

# evaluate accuracy
echo "evaluating accuracy..."
./accuracy-test --bam groot.bam --numReads $NUM_READS > ../accuracy-for-${NUM_READS}-reads.txt

# clean up
cd ..
rm -r ${NUM_READS}
echo "done"
