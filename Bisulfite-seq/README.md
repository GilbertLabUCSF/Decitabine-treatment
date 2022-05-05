# Decitabine DNA methylation data
I'm using previously published Bisulfite-Seq data. 
- [HL60 10ng DMSO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4518676)
- [HL60 10ng Decitabine](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4518677)

Download SRA fastq https://erilu.github.io/python-fastq-downloader/

## Preprocessing 

<!-- > MethPipe: a computational pipeline for analyzing bisulfite sequencing data ([link](http://smithlabresearch.org/software/methpipe/))

> Update July 2021: MethPipe now accepts SAM input after the read-mapping phase. Our old mr format is no longer supported. For short-read bisulfite mapping, we have a new tool called [abismal](https://github.com/smithlabcode/abismal/).

Okay! I'm using abismal. I installed `v1.0.0` through conda - https://anaconda.org/bioconda/abismal. 
 -->
__[(I) Bismark Genome Preparation](https://github.com/FelixKrueger/Bismark/tree/master/Docs#i-bismark-genome-preparation)__
```
bismark_genome_preparation genomes/hg38/gencode.v34/
```