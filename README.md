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
___
List of other repositories developed as part of this project:
- Integrated methods for meRIP-seq data analysis – https://github.com/abearab/imRIP
- Pathway analysis using PAGE algorithm and downstream analysis – https://github.com/abearab/pager
- Scripts for mapping NGS reads to HERVs – https://github.com/abearab/HERVs

___
This repo has been developed in Gilbert and Goodarzi labs at UCSF (2019-2022) by [Abolfazl Arab](https://github.com/abearab) and incldue some original/revised codes from  Alex Y. Ge,  Hosseinali Asgharian, and Hani Goodarzi.
