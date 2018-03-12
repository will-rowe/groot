<div align="center">
    <img src="paper/img/misc/groot-logo-with-text.png?raw=true?" alt="groot-logo" width="250">
    <h3><a style="color:#D5672C">G</a>raphing <a style="color:#D5672C">R</a>esistance <a style="color:#D5672C">O</a>ut <a style="color:#D5672C">O</a>f me<a style="color:#D5672C">T</a>agenomes</h3>
    <hr>
    <a href="https://travis-ci.org/will-rowe/groot"><img src="https://travis-ci.org/will-rowe/groot.svg?branch=master" alt="travis"></a>
    <a href='http://groot-documentation.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/groot-documentation/badge/?version=latest' alt='Documentation Status' /></a>
    <a href="https://goreportcard.com/report/github.com/will-rowe/groot"><img src="https://goreportcard.com/badge/github.com/will-rowe/groot" alt="reportcard"></a>
    <a href="https://github.com/will-rowe/groot/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
    <a href="https://zenodo.org/badge/latestdoi/117543539"><img src="https://zenodo.org/badge/117543539.svg" alt="DOI"></a>
    <a href="https://gitter.im/groot-help/Lobby"><img src="https://img.shields.io/badge/chat-gitter-8A2BE2.svg" alt="bioconda"></a>
    <a href="https://bioconda.github.io/recipes/groot/README.html"><img src="https://anaconda.org/bioconda/groot/badges/downloads.svg" alt="bioconda"></a>
</div>

***

## Overview

`GROOT` is a tool to type **Antibiotic Resistance Genes** (**ARGs**) in metagenomic samples. It combines variation graph representation of gene sets with an LSH indexing scheme to allow for fast classification of metagenomic reads. Subsequent hierarchical local alignment of classified reads against graph traversals facilitates accurate reconstruction of full-length gene sequences using a simple scoring scheme.

`GROOT` will output an ARG alignment file (in BAM format) that contains the graph traversals possible for each query read; the alignment file is then used to generate a simple ARG typing report.


## Installation

Check out the [releases](https://github.com/will-rowe/groot/releases) to download a binary. Alternatively, install using Bioconda or compile the software from source.

### Bioconda

```
conda install groot
```

> note: if using Conda make sure you have added the [Bioconda](https://bioconda.github.io/) channel first

### Source

`GROOT` is written in Go (v1.9) - to compile from source you will first need the [Go tool chain](https://golang.org/doc/install). Once you have it, try something like this to compile:

```bash
# Clone this repository
git clone https://github.com/will-rowe/groot.git

# Go into the repository and get the package dependencies
cd groot
go get -d -t -v ./...

# Compile the program
go build ./

# Call the program
./groot --help
```


## Quick Start

`GROOT` is called by typing **groot**, followed by the subcommand you wish to run. There are three main subcommands: **index**, **align** and **report**. This quick start will show you how to get things running but it is recommended to follow the [documentation](http://groot-documentation.readthedocs.io/en/latest/?badge=latest).

```bash
# Get a pre-clustered ARG database
groot get -d arg-annot

# Create graphs and index
groot index -i arg-annot.90 -o groot-index -l 100

# Align reads and report
groot align -i groot-index -f reads.fq | groot report
```
>note: you must run the index subcommand using your query read length (e.g. `-l 100` for 100bp reads)


## Further Information & Citing

Please [readthedocs](http://groot-documentation.readthedocs.io/en/latest/?badge=latest) for more extensive documentation and a [tutorial](https://groot-documentation.readthedocs.io/en/latest/tutorial.html).

We also have a preprint paper describing `GROOT`:

> [Rowe WPM, Winn MD. Indexed variation graphs for efficient and accurate resistome profiling. bioRxiv. 2018. doi: https://doi.org/10.1101/270835](https://doi.org/10.1101/270835)
