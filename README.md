<div align="center">
  <img src="paper/img/misc/groot-logo.png?raw=true?" alt="groot-logo" width="250">
  <br><br>
  <h1>GROOT</h1>
  <h3><a style="color:#FF6600">G</a>raphing <a style="color:#FF6600">R</a>esistance <a style="color:#FF6600">O</a>ut <a style="color:#FF6600">O</a>f me<a style="color:#FF6600">T</a>agenomes</h3>
  <hr>
  <a href="https://travis-ci.org/will-rowe/groot"><img src="https://travis-ci.org/will-rowe/groot.svg?branch=master" alt="travis"></a>
  <a href='http://groot-documentation.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/groot-documentation/badge/?version=latest' alt='Documentation Status' /></a>
 <a href="https://goreportcard.com/report/github.com/will-rowe/groot"><img src="https://goreportcard.com/badge/github.com/will-rowe/groot" alt="reportcard"></a>
  <a href="https://github.com/will-rowe/groot/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
  <a href="https://bioconda.github.io/recipes/groot/README.html"><img src="https://anaconda.org/bioconda/groot/badges/version.svg" alt="bioconda"></a>
  <a href="https://zenodo.org/badge/latestdoi/117543539"><img src="https://zenodo.org/badge/117543539.svg" alt="DOI"></a>
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

> If using Conda, make sure you have added the [Bioconda](https://bioconda.github.io/) channel first.

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

`GROOT` is called by typing **groot**, followed by the command you wish to run. There are three main commands: **index**, **align** and **report**. This quick start will show you how to get things running but it is recommended to follow the [documentation](http://groot-documentation.readthedocs.io/en/latest/?badge=latest).

```
# get a pre-clustered ARG database
groot get -d arg-annot

# create graphs and index
groot index -i arg-annot.90 -o groot-index

# align reads and report
groot align -i groot-index -f reads.fq | groot report
```


## Further Information and citing

Please [readthedocs](http://groot-documentation.readthedocs.io/en/latest/?badge=latest) for more extensive documentation and a [tutorial](https://groot-documentation.readthedocs.io/en/latest/tutorial.html).

The paper describing `GROOT` is in preparation and will be pre-printed soon.
