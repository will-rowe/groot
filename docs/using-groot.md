# Using GROOT

**GROOT** uses a subcommand syntax; call the main program with `groot` and follow it by the subcommand to indicate what  action to take.

This page will cover a worked example, details on the available GROOT commands and some tips for using the program.

For more information on the graphs that **GROOT** uses, please read the [groot-graphs](https://groot-documentation.readthedocs.io/en/latest/groot-graphs.html) page.

***

## An example

This example will take us from raw reads to resistome profile for a single metagenome

Get some sequence data:
```
fastq-dump SRR4454613
```

Get a pre-clustered ARG database:
```
groot get -d resfinder
```

Create variation graphs from the ARG-database and index:
```
groot index -i resfinder.90 -o groot-index -p 8 -l 100
```

Align the reads and report ARGs
```
groot align -i groot-index -f SRR4454613.fastq -p 8 | groot report -c 1
```

***

## GROOT subcommands

### get

The ``get`` subcommand is used to download a pre-clustered ARG database that is ready to be indexed. Here is an example:

```
groot get -d resfinder -i 90 -o .
```

The above command will download the [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) database, which has been clustered at 90% identity, check its integrity and unpacks it to the current directory. The database will be named ``resfinder.90`` and contain several ``*.msa`` files.

Flags explained:

* ``-d``: which database to get
* ``-i``: the identity threshold at which the database was clustered (note: only 90 currently available)
* ``-o``: directory to save the database to

The following databases are available:

* arg-annot (default)
* resfinder
* card
* groot-db
* groot-core-db

These databases were clustered by sequence identity and stored as **Multiple Sequence Alignments** (MSAs). See [groot-databases](https://groot-documentation.readthedocs.io/en/latest/groot-databases.html) for more info.

### index

The ``index`` subcommand is used to create variation graphs from a pre-clustered ARG database and then index them. Here is an example:

```
groot index -i resfinder.90 -o groot-index -p 8 -l 100
```

The above command will create a variation graph for each cluster in the ``resfinder.90`` database and initialise an LSH forest index. It will then move through each graph traversal using a 100 node window, creating a MinHash signature for each window. Finally, each signature is added to the LSH Forest index. The index will be named ``groot-index`` and contain several files.

Flags explained:

* ``-i``: which database to use
* ``-o``: where to save the index
* ``-p``: how many processors to use
* ``-l``: length of window to use (**should be similar to the length of query reads (+/- 10 bases)**)

The same index can be used on multiple samples, however, it should be re-indexed if you are using different read lengths or wish to change the seeding parameters.

Some more flags that can be used:

* ``-j``: Jaccard similarity threshold for seeding a read (used to calibrate LSH Forest index)
* ``-k``: size of k-mer to use for MinHashing
* ``-s``: length of MinHash signature

> Important: GROOT can only accept MSAs as input. You can cluster your own database or use `groot get` to obtain a pre-clustered one.

### align

The ``align`` subcommand is used to align reads against the indexed variation graphs. Here is an example:

```
groot align -i groot-index -f file.fastq -p 8 > ARG-reads.bam
```

The above command will seed the fastq reads against the indexed variation graphs. It will then perform a hierarchical local alignment of each seed against the variation graph traversals. The output alignment is essentially the ARG classified reads (which may be useful) and can then be used to report full-length ARGs (using the `report` subcommand).

Flags explained:

* ``-i``: which index to use
* ``-f``: what FASTQ reads to align
* ``-p``: how many processors to use

Data can streamed in and out of the align subcommand. For example:

```
gunzip -c file.gz | ./groot align -i groot-index -p 8 | ./groot report
```
> this feature is being worked on and will bring more utility in future releases.

Multiple FASTQ files can be specified as input, however all are treated as the same sample and paired-end info isn't used. To specify multiple files, make sure they are comma separated (``-i fileA.fq,fileB.fq``).

Some more flags that can be used:

* ``--trim``: enable quality based trimming of read prior to alignment
* ``-q``: minimum base quality during trimming
* ``-l``: minimum read length post trimming
* ``-c``: maximum number of bases to clip from 3' end during local alignment
* ``-o``: directory to save variation graphs to

### report

The ``report`` subcommand is used to processes graph traversals and generate a resistome profile for a sample. Here is an example:

```
groot report -i ARG-reads.bam -c 1
```

Flags explained:

* ``-i``: the input BAM file (output from `groot align` subcommand)
* ``-c``: the coverage needed to report an ARG (e.g. 0.95 = 95% ARG bases covered by reads)

Some more flags that can be used:

* ``--lowCov``: overrides `c` option and will report ARGs which may not be covered at the 5'/3' ends
* ``--plotCov``: outputs coverage plot for each ARG reported

***

## Tips

Will update soon...
