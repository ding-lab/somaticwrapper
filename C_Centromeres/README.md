# Create BED files defining centromeres and other featues

A centromere BED file is used for Pindel to define regions to be excluded 
from analysis.  While optional, it is used for performance purposes.

Three such files are distributed:
* pindel-centromere-exclude.bed  -  this file is for GRCh37, and was used for internal development purposes.  It has an unknown provenance.
* ucsc-centromere.GRCh37.bed - Data generated from UCSC Table Browser for GRCh37 assembly as described below
* ucsc-centromere.GRCh38.bed - Data generated from UCSC Table Browser for GRCh38 assembly as described below

There is no reason in most cases to do anything in this directory, but instructions of how to manually obtain UCSC data are below.

## ucsc-centromere.GRCh37.bed

Table is obtained manually from [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables)
* *Assembly:* GRCh37
* *Group:* All Tables
* *Table:* Gap
* *Filter:* Type matches "centromere"
* *Output format:* BED
* *Output File:* ucsc-centromere.GRCh37.bed
* Defaults for *Output gap as BED*


## ucsc-centromere.GRCh38.bed

Table is obtained manually from [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables).  Note that GRCh38 treats centromeres differently than GRCh37, as described [here](https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/SaR2y4UNrWg).
* *Assembly:* GRCh38
* *Group:* All Tables
* *Table:* Centromeres
* *Output format:* BED
* *Output File:* ucsc-centromere.GRCh38.bed
* Defaults for *Output gap as BED*
