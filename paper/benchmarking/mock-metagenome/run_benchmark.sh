#!/bin/env bash
#
# This runs the mock metagenome benchmark
#
# The CAMI metagenomic dataset is spiked with ARGs from the CARD database
# Randomreads are generated at the specified coverage level (errors are allowed)
# Groot and ARGOAP are both run
#
# run with parallel using a command like:
# cat test-coverage-values | parallel --gnu "sh run_benchmark.sh {}"
#
# Requires:
#
# wget
# bioawk (conda install bioawk)
# BBMap
# Go 1.9
# gtime (brew install gnu-time)

# test parameters
COVERAGE=$1
THREADS=1
READ_LEN=150
K_SIZE=7
SIG_SIZE=128
JS=0.99
ARG_SAMPLE_SIZE=10
echo "test parameters:"
echo "  - coverage: " $COVERAGE
echo "  - read length: " $READ_LEN
echo "  - k-mer size: " $K_SIZE
echo "  - signature size: " $SIG_SIZE
echo "  - Jaccard similarity: " $JS

mkdir coverage-${COVERAGE} && cd coverage-${COVERAGE}

# download CAMI genomes (low complexity)
echo "downloading CAMI genomes..."
wget https://storage.googleapis.com/cami-data-eu/CAMI_low/source_genomes_low.tar.gz > /dev/null 2>&1
# wget https://storage.googleapis.com/cami-data-eu/CAMI_medium/source_genomes_medium.tar.gz > /dev/null 2>&1
# wget https://storage.googleapis.com/cami-data-eu/CAMI_high/source_genomes_high.tar.gz > /dev/null 2>&1
tar -xvf source_genomes_low.tar.gz  > /dev/null 2>&1 && cd source_genomes
for i in *.f*; do
fuse.sh in=$i out=$i.fused.fna pad=0 > /dev/null 2>&1
rm $i
done

# randomly sample some ARGs from the CARD database
echo "sampling CARD genes..."
bioawk -c fastx -v k=${ARG_SAMPLE_SIZE} '{y=x++<k?x-1:int(rand()*x);if(y<k)a[y]=">"$name"\n"$seq}END{for(z in a)print a[z]}' ../../../../data/full-ARG-databases/card/card-1.1.2.fna > randomARGs.fna

# create artificial metagenome
echo "creating artificial metagenome..."
tail -n +1 *.fna > genomes.fna
sed '/^==>/ d' < genomes.fna > ../groot-reference-data.fna
cd ..
randomreads.sh ref=groot-reference-data.fna len=${READ_LEN} out=groot-reads.fq metagenome=true coverage=${COVERAGE} adderrors=true maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 > /dev/null 2>&1
split -l $(($(wc -l < groot-reads.fq) / 2)) groot-reads.fq split
mv splitaa groot-reads_1.fq
mv splitab groot-reads_2.fq
rm groot-reads.fq

# install GROOT
echo "installing GROOT..."
go build ../../../..

# index the CARD database
echo "indexing the CARD database..."
gtime -f "\tmax. resident set size (kb): %M\n\tCPU usage: %P\n\ttime (wall clock): %E\n\ttime (CPU seconds): %S\n" ./groot index -i ../../../data/clustered-ARG-databases/card-90 -o groot-index -l ${READ_LEN} -k ${K_SIZE} -s ${SIG_SIZE} -j ${JS} -p ${THREADS}

# run the alignment and generate report
echo "aligning the test reads..."
gtime -f "\tmax. resident set size (kb): %M\n\tCPU usage: %P\n\ttime (wall clock): %E\n\ttime (CPU seconds): %S\n" ./groot align -i groot-index -f groot-reads_1.fq,groot-reads_2.fq -p $THREADS > groot-out.bam
echo "running groot report..."
./groot report -i groot-out.bam -c 1 -p $THREADS > groot.report

# process the report and check against input ARGs
cut -f 1 groot.report > groot-foundARGs.txt
found=0
fp=0
while read line; do
if grep -q "${line}" groot-reference-data.fna
then
found=$((found + 1))
else
fp=$((fp + 1))
fi
done <groot-foundARGs.txt
echo "GROOT finished - results are in 'results-for-groot.txt'"
echo "GROOT results:" > results-for-groot.txt
echo "correctly found ${found} / ${ARG_SAMPLE_SIZE} ARGs" >> results-for-groot.txt
echo "also found ${fp} false positives" >> results-for-groot.txt

# run ARGsOAP (stage 1)
echo "running ARGsOAP (stage 1)..."
git clone  https://github.com/biofuture/Ublastx_stageone.git > /dev/null 2>&1
cd Ublastx_stageone
printf "SampleID\tName\tCategory\n1\tgroot-reads\tblank\n" > input_files.txt
./ublastx_stage_one -i ../ -o argsOAP-output -m input_files.txt > /dev/null 2>&1
mkdir ../argOAP-upload
mv argsOAP-output/extracted.fa ../argOAP-upload/ && mv argsOAP-output/meta_data_online.txt ../argOAP-upload/
cd ..
echo "argsOAP stage 1 finished - upload the 2 files in the 'argOAP-upload' directory to http://smile.hku.hk/SARGs and run argsOAP stage 2"

# clean up
yes | rm -r groo* source_genome* ref Ublastx_stageone
echo "done"
