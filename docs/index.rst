Welcome to GROOT's wiki!
============================

`GROOT` is a tool to profile **Antibiotic Resistance Genes** (**ARGs**) in metagenomic samples.

The main advantages of `GROOT` over existing tools are:

* quick classification of reads to candidate ARGs
* accurate annotation of full-length ARGs
* can run on a laptop in minutes

`GROOT` aligns reads to ARG variation graphs, producing an alignment file that contains the graph traversals possible for each query read. The alignment file is then used to generate a simple resistome profile report.

----

Overview
------------------------------------------------------

Antimicrobial resistance remains a major threat to global health. Profiling the collective antimicrobial resistance genes within a metagenome (the "resistome") facilitates greater understanding of antimicrobial resistance gene diversity and dynamics. In turn, this can allow for gene surveillance, individualised treatment of bacterial infections and more sustainable use of antimicrobials. However, resistome profiling can be complicated by high similarity between reference genes, as well as the sheer volume of sequencing data and the complexity of analysis workflows. We have developed an efficient and accurate method for resistome profiling that addresses these complications and improves upon currently available tools.

Our method combines variation graph representation of gene sets with an LSH Forest indexing scheme to allow for fast classification of metagenomic reads using similarity-search queries. Subsequent hierarchical local alignment of classified reads against graph traversals facilitates accurate reconstruction of full-length gene sequences using a scoring scheme. GROOT runs on a laptop and can process a typical 2 gigabyte metagenome in 2 minutes using a single CPU.

----

Contents
------------------------------------------------------
.. toctree::
   :maxdepth: 2

   using-groot
   groot-graphs
   groot-databases
   tutorial
