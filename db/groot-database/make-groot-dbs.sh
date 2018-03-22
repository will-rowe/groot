#!/bin/env bash
#
# This script generates the groot-core-db and the groot-db
#
# It downloads the latest versions of the following databases:
#
#   - argannot
#   - resfinder
#   - card
#   - megares
#
# It then either:
#   1. groot-db:  merge the databases, remove duplicates and then cluster
#   or
#   2. groot-core-db: extract common ARGs from all databases and then cluster these common ARGs
#
# REQUIRES: Vsearch, SeqKit
#

echo "making the groot and groot-core databases..."
mkdir tmp && cd $_

# Download the latest CARD database
mkdir card && cd $_
wget --no-check-certificate -qO- https://card.mcmaster.ca/download | grep -Eo "download/0/broadstreet[a-zA-Z0-9./?=_-]*" | sort | uniq | tail -n 1 > card-db-version
cardLink=$(sed  's/^/https:\/\/card.mcmaster.ca\//g' card-db-version)
wget --no-check-certificate -O card-db.tar.gz $cardLink
tar -xvf card-db.*
awk '/>/{sub(">",">groot-db_CARD__")}1' nucleotide_fasta_protein_homolog_model.fasta > ../card-refs.fna
cd .. && rm -r card

# Download the latest ARG-annot database (V3)
wget http://en.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/1424/arg-annot-nt-v3-march2017_doc.fasta -O argannot-refs.fna
awk '/>/{sub(">",">groot-db_ARGANNOT__")}1' argannot-refs.fna > tmp && mv tmp argannot-refs.fna

# Download the latest ResFinder database
mkdir resfinder && cd $_
wget https://bitbucket.org/genomicepidemiology/resfinder_db/get/dc33e2f9ec2c.zip -O resfinder.zip
unzip resfinder.zip
awk 'FNR==1{print ""}1' genomic*/*.fsa > resfinder-refs.fna
awk '/>/{sub(">",">groot-db_RESFINDER__")}1' resfinder-refs.fna > ../resfinder-refs.fna
cd .. && rm -r resfinder

# Download the latest megres database
#mkdir megares && cd $_
#wget --no-check-certificate -qO- https://megares.meglab.org/download/index.php | grep -Eo "megares_v.*/megares_database.*[0-9].fasta" | sort | uniq | tail -n 1 > megres-db-version
#megaresLink=$(sed  's/^/https:\/\/megares.meglab.org\/download\//g' megres-db-version)
#wget --no-check-certificate -O ../megares-refs.fna $megaresLink
#cd .. && rm -r megares

# Create a reference file for the core database
seqkit common *.fna --by-seq --ignore-case -o core-args.fasta -j 8

# Cluster core set
mkdir groot-core-db.90 && cd $_
vsearch --cluster_size ../core-args.fasta --id 0.90 --msaout MSA.tmp
awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.tmp && rm MSA.tmp
date +%x_%H:%M:%S:%N | sed 's/\(:[0-9][0-9]\)[0-9]*$/\1/' > timestamp.txt
cd ..

# Create a reference file for the complete database
cat *.fna > all-args.fasta
seqkit rmdup --by-seq --ignore-case -j 8 -o all-args.dedup.fasta < all-args.fasta

# Cluster total set
mkdir groot-db.90 && cd $_
vsearch --cluster_size ../all-args.dedup.fasta --id 0.90 --msaout MSA.tmp
awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.tmp && rm MSA.tmp
date +%x_%H:%M:%S:%N | sed 's/\(:[0-9][0-9]\)[0-9]*$/\1/' > timestamp.txt
cd ..

# Finish up
mv groot* ..
cd .. && rm -r tmp
