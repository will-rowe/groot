#!/bin/bash
set -ex

# samtools 1.5
wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
tar -xvf samtools-1.5.tar.bz2
cd samtools-1.5 && ./configure && make
sudo make PREFIX=/usr install
cd ..

# bamUtil
wget --no-check-certificate https://github.com/statgen/libStatGen/archive/v1.0.14.tar.gz
tar -xvf v1.0.14.tar.gz
cd libStatGen-1.0.14
make all
cd ..
wget --no-check-certificate https://github.com/statgen/bamUtil/archive/v1.0.14.tar.gz -O bamUtil.tar.gz
tar -xvf bamUtil.tar.gz
cd bamUtil-1.0.14
make all LIB_PATH_BAM_UTIL=../libStatGen-1.0.14/ PREFIX=/usr
cd
