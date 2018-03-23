# GROOT Databases

This is a brief overview of the databases used by **GROOT**.

***

## Overview

As mentioned on the [groot-graphs](/groot-graphs.html) page, **GROOT** creates variation graphs from [Multiple Sequence Alignments](https://en.wikipedia.org/wiki/Multiple_sequence_alignment) (**MSAs**) and stores them as `groot graphs`.

The **MSAs** represent a clustered database, each cluster is a collection of sequences which share high nucleotide identity. **GROOT** reads in each **MSA** file and converts it to a graph.

A set of clustered databases is provided (use the `groot get` subcommand) or you can generate your own clustered database. To cluster a database yourself, you can use the following commands on any multifasta file containing the sequences you want to use with **GROOT**:

``` bash
mkdir CLUSTERED-DB && cd $_

vsearch --cluster_size /path/to/ARGs.fna --id 0.90 --msaout MSA.tmp

awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.tmp && rm MSA.tmp

cd ..
```

* the above snippet will create a clustered database in the directory CLUSTERED-DB - now you can run `groot index`

```
groot index -i ./CLUSTERED-DB
```


## groot-db and groot-core-db

As mentioned earlier, the `groot get` subcommand can download a pre-clustered database for you to use with **GROOT**. The following databases are available:

* arg-annot (default)
* resfinder
* card
* groot-db
* groot-core-db

The `groot-db` and `groot-core-db` are both databases that are derived from ResFinder, ARG-annot and CARD. They have been included after requests by several users for a combination of available **ARG** databases. They were made as follows:

`groot-db` is made by combining all sequences in ResFinder, ARG-annot and CARD. Duplicates are removed and the sequences are then clustered at 90% identity.

`groot-core-db` is made by combining sequences that are present in each of the ResFinder, ARG-annot and CARD databases. One copy of each sequence is kept and then this collection is clustered at 90% identity.

Both `groot-db` and `groot-core-db` prepend a tag to each reference sequence so that the origin can be determined. For example:

```
>*groot-db_ARGANNOT__(AGly)APH-Stph:HE579073:1778413-1779213:801
```

In the directory downloaded by `groot get`, there will also be a timestamp that tells you when the database was created. The **GROOT** database will have used the most recently available versions of ResFinder/CARD/ARG-annot. The script used to do this is available in the **GROOT** repo:

```
db/groot-database/make-groot-dbs.sh
```
