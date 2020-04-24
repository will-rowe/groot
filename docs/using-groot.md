# Using GROOT

**GROOT** uses a subcommand syntax; call the main program with `groot` and follow it by the subcommand to indicate what action to take.

This page will cover a worked example, details on the available GROOT commands and some tips for using the program.

For more information on the graphs that **GROOT** uses, please read the [groot-graphs](https://groot-documentation.readthedocs.io/en/latest/groot-graphs.html) page.

---

## An example

This example will take us from raw reads to resistome profile for a single metagenome

Get some sequence data:

```sh
fastq-dump SRR4454613
```

Get a pre-clustered ARG database:

```sh
groot get -d resfinder
```

Create variation graphs from the ARG-database and index:

```sh
groot index -m resfinder.90 -i grootIndex -w 100 -p 8
```

Align the reads and report ARGs

```sh
groot align -i grootIndex -f SRR4454613.fastq -p 8 | groot report -c 0.95
```

---

## GROOT subcommands

### get

The `get` subcommand is used to download a pre-clustered ARG database that is ready to be indexed. Here is an example:

```
groot get -d resfinder -o . --identity 90
```

The above command will download the [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) database, which has been clustered at 90% identity, check its integrity and unpacks it to the current directory. The database will be named `resfinder.90` and contain several `*.msa` files.

Flags explained:

- `-d`: which database to get
- `-o`: directory to save the database to
- `--identity`: the identity threshold at which the database was clustered (note: only 90 currently available)

The following databases are available:

- arg-annot (default)
- resfinder
- card
- groot-db
- groot-core-db

These databases were clustered by sequence identity and stored as **Multiple Sequence Alignments** (MSAs). See [groot-databases](https://groot-documentation.readthedocs.io/en/latest/groot-databases.html) for more info.

### index

The `index` subcommand is used to create variation graphs from a pre-clustered ARG database and then index them. Here is an example:

```
groot index -m resfinder.90 -i grootIndex -w 100 -p 8
```

The above command will create a variation graph for each cluster in the `resfinder.90` database and initialise an LSH forest index. It will then move through each graph traversal using a 100 node window, creating a MinHash sketch for each window. Finally, each sketch is added to the LSH Forest index. The index will be named `grootIndex` and contain several files.

Flags explained:

- `-m`: a directory of MSA files (the database from `groot get`)
- `-i`: where to save the index
- `-w`: the window size to use (**should be similar to the length of query reads**)
- `-p`: how many processors to use

The same index can be used on multiple samples, however, it should be re-indexed if you wish to change the seeding parameters.

Some more flags that can be used:

- `-k`: size of k-mer to use for MinHashing
- `-s`: size of MinHash sketch
- `-x`: number of partitions in the LSH Ensemble index
- `-y`: maxK in the LSH Ensemble index

> Important: GROOT can only accept MSAs as input. You can cluster your own database or use `groot get` to obtain a pre-clustered one.

### align

The `align` subcommand is used to align reads against the indexed variation graphs. Here is an example:

```
groot align -i grootIndex -f file.fastq -t 0.97 -p 8 > ARG-reads.bam
```

The above command will seed the fastq reads against the indexed variation graphs. It will then perform a hierarchical local alignment of each seed against the variation graph traversals. The output alignment is essentially the ARG classified reads (which may be useful) and can then be used to report full-length ARGs (using the `report` subcommand).

Flags explained:

- `-i`: which index to use
- `-f`: what FASTQ reads to align
- `-t`: the containment threshold for seeding reads
- `-p`: how many processors to use

Data can streamed in and out of the align subcommand. For example:

```
gunzip -c file.gz | ./groot align -i grootIndex -p 8 | ./groot report
```

Multiple FASTQ files can be specified as input, however all are treated as the same sample and paired-end info isn't used. To specify multiple files, make sure they are comma separated (`-f fileA.fq,fileB.fq`) or use gunzip/cat with a wildcard (gunzip -c \*.fq.gz | groot...).

### report

The `report` subcommand is used to processes graph traversals and generate a resistome profile for a sample. Here is an example:

```
groot report --bamFile ARG-reads.bam -c 1
```

Flags explained:

- `--bamFile`: the input BAM file (output from `groot align` subcommand)
- `-c`: the coverage needed to report an ARG (e.g. 0.95 = 95% ARG bases covered by reads)

Some more flags that can be used:

- `--lowCov`: overrides `c` option and will report ARGs which may not be covered at the 5'/3' ends
