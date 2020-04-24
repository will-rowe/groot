#!/bin/env bash
#
# This test checks that the correct gene is reported by the GROOT workflow.
#
set -e

# test parameters
TESTDIR=tmp-for-groot-travis
READS=../data/bla-b7-150bp-5x.fq
THREADS=1
READ_LEN=150
K_SIZE=31
SIG_SIZE=42
CT=0.99
COV=0.97

# make a test dir
mkdir $TESTDIR && cd $TESTDIR

# build the prog
go build -o groot ../../

# get the db
echo "downloading the ARG-annot database..."
./groot get -d arg-annot > /dev/null 2>&1

# index the ARGannot database
echo "indexing the ARG-annot database..."
#gtime -f "\tmax. resident set size (kb): %M\n\tCPU usage: %P\n\ttime (wall clock): %E\n\ttime (CPU seconds): %S\n" \
./groot index -m arg-annot.90 -i index -w $READ_LEN -k $K_SIZE -s $SIG_SIZE -p $THREADS

# align the test reads
echo "aligning reads..."
#gtime -f "\tmax. resident set size (kb): %M\n\tCPU usage: %P\n\ttime (wall clock): %E\n\ttime (CPU seconds): %S\n" \
./groot align -i index -f $READS -p $THREADS -t $CT > groot.bam

# report
echo "reporting..."
./groot report --bamFile groot.bam -c $COV > groot.report

# check that bla-b7 is the only arg reported
echo "checking..."
numReportedARGs=`wc groot.report | awk '{print $1}'`
if [[ $numReportedARGs == "0" ]]; then
    echo "failed: no ARGs reported by GROOT";
    exit 1; 
fi
if [[ $numReportedARGs != "1" ]]; then
    echo "failed: too many ARGs reported by GROOT - " $numReportedARGs;
    exit 1; 
fi
reportedARG=`cut -f 1 groot.report`
if [[ $reportedARG != "argannot~~~(Bla)B-7~~~AF189304:1-747" ]]; then
    echo "failed: GROOT got the wrong gene - expected (Bla)B-7 but got " $reportedARG
    exit 1;
fi

# clean up
cd ..
rm -r $TESTDIR
echo "passed"
