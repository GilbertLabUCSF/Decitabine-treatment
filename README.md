**A multiomics approach reveals RNA dynamics promote cellular sensitivity to DNA hypomethylation**

This repository contains all the scripts necessary for producing the results of the manuscript, for which a preprint is available at:

[![DOI:10.1101/2022.12.14.518457](http://img.shields.io/badge/DOI-10.1101/2022.12.14.518457-B31B1B.svg)](https://www.biorxiv.org/content/early/2022/12/14/2022.12.14.518457)

Raw data is also available at:

[![GEO:GSE222886](https://img.shields.io/badge/GEO-GSE222886-green.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222886)

___
List of notebooks and scripts used for the analysis in the manuscript is as follows:
```bibtex
├── DAC
│   ├── Bisulfite-seq
│   │   ├── 0-preprocessing.ipynb
│   │   ├── 1-differential-DNAme-analysis.ipynb
│   ├── CRISPRi-screen
│   │   ├── CRISPRi-screen-cell-line-consistency.ipynb
│   │   ├── CRISPRi-screen-geneset-analysis-gsea.ipynb
│   │   ├── CRISPRi-screen-volcano-plot.ipynb
│   ├── RNA-seq
│   │   ├── exp-differential-expression.ipynb
│   │   ├── pathway-enrichment.ipynb
│   │   ├── stbl-differential-stability.ipynb
│   │   ├── stbl-run-REMBRANDTS.ipynb
│   ├── Ribo-seq
│   │   ├── ribolog_delta_te.ipynb
│   │   └── te-downstream.ipynb
│   └── meRIP-seq
│       ├── coverage-plots-herv.ipynb
│       ├── coverage-plots.ipynb
│       ├── exomepeak-call-peaks.ipynb
│       ├── radar-diff-peaks-herv.ipynb
│       ├── radar-diff-peaks.ipynb
│       ├── radar-volcano-plot.ipynb
├── DAC-rg3039
│   ├── RNA-seq
│   │   ├── ERV-differential-expression.ipynb
│   │   ├── ERV-upsetplot.ipynb
│   │   ├── drug-comb-differential-expression.ipynb
│   │   ├── drug-comb-downstream-plots.ipynb
│   │   ├── drug-comb-geneset-enrichment.ipynb
│   │   ├── drug-comb-pathway-enrichment.ipynb
│   │   ├── drug-comb-scatter-plot.ipynb
│   │   └── scripts
│   │       ├── deseq_drug_combo_pipeline.R
│   │       └── salmon.sh
│   ├── apoptosis_cell_cycle_plots.ipynb
│   └── synergy.ipynb
├── combined_analysis
│   ├── cell-line-consistency.ipynb
│   ├── intersectional-analysis.ipynb
└── scripts
    ├── fastqc.sh
    ├── util.R
    └── util.py

```
Pipelines

Differential gene expression analysis: The Salmon-tximport-DESeq2 pipeline
> We used a workflow hereafter referred to as the “Salmon-tximport-DESeq2 pipeline” to perform differential gene expression analysis. Salmon (version 1.2.1) was first used to quantify transcript abundance53. A Salmon index was generated using the GENCODE75 (version 34) genome annotation, and subsequently the `salmon quant` tool was used with the `--validateMappings` option to calculate transcript abundances. Then, the R package tximport52 was used to import Salmon results into R and perform data preparation. The `summarizeToGene` function was used to collapse transcript abundances to the gene level. From here, the R package DESeq251 was used for differential gene expression analysis. We first extracted normalized counts for each RNA-seq experiment using DESeq2 by running the `estimateSizeFactors` function and then the `counts` function with option `normalized=TRUE`. For each individual experiment, the DESeq2 statistical model was modified based on the experimental design. For experimental designs with multiple variables (e.g., multiple drug conditions, time points, etc.), we used the likelihood ratio test (LRT) to perform differential expression analysis. The LRT is conceptually similar to an analysis of variance (ANOVA) calculation in a linear regression model76. In these cases, we specified the model design in the `DESeq2` function as `~0 + variable1 + variable2 + variable1:variable2` and the argument `test=LRT`. In simple experimental designs with one variable (e.g., DMSO vs. decitabine treatment), DESeq2 was used with default options (i.e., a Wald test was used instead of a LRT). In these cases, the model design was specified as `~cond`. For experiments with batch effects, the model design was specified as `~cond + reps`.

Differential RNA stability analysis: The STAR-featureCounts-REMBRANDTS-limma pipeline
> For analyses which required measurements of pre-mRNA and mature mRNA abundances from RNA-seq samples (i.e., differential RNA stability analysis), we used a workflow hereafter referred to as the “STAR-featureCounts-REMBRANDTS-limma pipeline”. RNA-seq sequencing reads were first aligned to the hg38 reference genome using STAR54 (version 2.7.3a). Then, featureCounts77 was used to quantify intron and exon level counts. Finally, REMBRANDTS was used to calculate mRNA stability as previously described56 (https://github.com/csglab/REMBRANDTS). Briefly, the package estimates a gene-specific bias function that is subtracted from Δexon–Δintron calculations to provide unbiased mRNA stability measurements. To assess differential RNA stability changes, we used limma55, which was designed for microarray experiments and serves a similar function to DESeq2, though it supports negative values (relevant for RNA stability analysis). The model designs used here are analogous to the designs for differential expression analysis described above.

___
List of other repositories developed as part of this project:
- Integrated methods for meRIP-seq data analysis – https://github.com/abearab/imRIP
- Pathway analysis using PAGE algorithm and downstream analysis – https://github.com/abearab/pager
- Scripts for mapping NGS reads to HERVs – https://github.com/abearab/HERVs

___
This repo has been developed in Gilbert and Goodarzi labs at UCSF (2019-2022) by [Abolfazl Arab](https://github.com/abearab) and incldue some original/revised codes from  Alex Y. Ge,  Hosseinali Asgharian, and Hani Goodarzi.
