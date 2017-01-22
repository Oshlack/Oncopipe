# Oncopipe

## Synopsis
Clinical pipeline for detecting gene fusions in RNAseq data utilising JAFFA.  The pipeline optionally calls a classifier for B-Cell Acute Lymphocytic Leukaemia (based on AllSorts).  As well as variant calling pipeline using Picard,GATK and VEP.
### Prerequisites
```
[AllSorts](https://github.com/Oshlack/AllSorts)
bbmap
blat
bowtie2
bpipe
fastqc
GATK
[JAFFA](https://github.com/Oshlack/JAFFA)
java
Picard
python
samtools v1.1
trimmomatic
RNA-Seqc
VEP

```
## Code Example
On meerkat:
```
module load java
module load bpipe
./pipeline/scripts/create_batch.sh {BATCH} B_ALL designs/B_ALL/B_ALL.bed

cd ./batches/{BATCH}/analysis
bpipe run ../../../pipeline.groovy ../samples.txt

```

## Motivation

## Installation

## Authors
* **Rebecca Louise Evans** - *Initial work* - [beccyl](https://github.com/beccyl)

## License
This project is licensed under GNU General Public License - see the [LICENSE.md](License.md) file for details.

## Acknowledgments
* [ssadedin](https://github.com/ssadedin) - provided the framework for bpipe and cpipe
* [nadiadavidson](https://github.com/nadiadavidson) - provided the underlying JAFFA "Just another fusion finding algorithm" on which this project is based.
* [Quarkins](https://github.com/Quarkins) - provided the algorithm for random forest classifier
* [JovMaksimovic](https://github.com/JovMaksimovic) - generic bpipe pipeline for rna-seq
* [aliciaoshlack](https://github.com/aliciaoshlack) - Team Leader
