here I'm using a custom [scallop](https://github.com/Kingsford-Group/scallop) + [gtfmerge](https://github.com/Kingsford-Group/rnaseqtools#gtfmerge) pipeline to create merged gtf annotation of hg19 and tinat transcripts introduced in this 

Alex
> I was curious myself so I decided to align our RNA-seq data to their chimeric transcripts... it looks promising, but let me know if there's anything in my analysis I need to change. Still working on seeing if we can identify alternative TSS without CAGE-seq. 
I'm guessing we can look for novel transcripts using [StringTie](https://ccb.jhu.edu/software/stringtie/) or cufflinks?

Ray
> I'm asking because in [this paper](https://www.nature.com/articles/ng.3889) they treated lung cancer cells with decitabine and performed Cage-seq to map and identify cryptic transcription start sites. They found cryptic TSSs that generated spliced variants, half of which were in-frame isoforms (and some of which resulted in chimeric transcripts, others original or truncated). This got me thinking about our project and whether we could see these chimeric/fusion transcripts (that are skewed toward the 5'end due Cage-seq and being TSS-induced) in HL60s.  


Luke
> Hi Ray, Alex and Abe
This is an interesting paper which claims genetic alterations to DNA methylation core genes drive dysregulated hematopoietic development due to CpG changes that regulate transcription factor binding (figure 4 and 5). How many TFs are methylation sensitive has been controversial in my understanding but is an interesting idea with lots of data within mylomonocytic cells. https://www.nature.com/articles/s41588-020-0595-4
I am not suggesting a specific plan but this is an interesting idea.

## HL-60 time-series Decitabine treatment RNA-Seq experiments
In order to test for any differences over multiple time points, once can use a design including the time factor, and then test using the **likelihood ratio test (LRT)**. Here, as we have control (DMSO) and treatment (Decitabine) time series, design formula containing the condition factor, the time factor, and the interaction of the two. In this case, using the likelihood ratio test with a reduced model which does not contain the interaction terms will test whether the condition induces a change in gene expression at any time point after the reference level time point (time 0).
(see [DESeq2 Time-series-experiments](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments) for more details)
<table>
  <tr>
    <td><img src=Volcano_plot.png width='500'>
    Volcano plot
  <tr>
<table>

## Heatmap clustering
<img src=Heatmap_clustering.png width='500'>

