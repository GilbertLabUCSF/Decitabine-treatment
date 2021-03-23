## HL-60 time-series Decitabine treatment RNA-Seq experiments
In order to test for any differences over multiple time points, once can use a design including the time factor, and then test using the **likelihood ratio test (LRT)**. Here, as we have control (DMSO) and treatment (Decitabine) time series, design formula containing the condition factor, the time factor, and the interaction of the two. In this case, using the likelihood ratio test with a reduced model which does not contain the interaction terms will test whether the condition induces a change in gene expression at any time point after the reference level time point (time 0).
(see [DESeq2 Time-series-experiments](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments) for more details)
<table>
  <tr>
    <td><img src=plots/PCA_filtered.png width='500'>
    PCA
    <td><img src=plots/Volcano_plot.png width='500'>
    Volcano plot
  <tr>
<table>

## Heatmap clustering
<img src=plots/Heatmap_clustering.png width='500'>

# Enrichment analysis
<table>
  <tr>
  <h2>human_ensembl<h2>
    <td><img src=6h_delta_exp/human_ensembl.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl.R.png width='500'>
    120h_delta_exp
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_encode_tf<h2>
    <td><img src=6h_delta_exp/human_ensembl_encode_tf.all.png width='500'>
    <td><img src=72h_delta_exp/human_ensembl_encode_tf.all.png width='500'>
    <td><img src=120h_delta_exp/human_ensembl_encode_tf.all.png width='500'>
  <tr>
<table> -->
<!-- <table>
  <tr>
  <h2>human_ensembl_msigdb_c1<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c1.all.png width='500'>
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c1.all.png width='500'>
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c1.all.png width='500'>
  <tr>
<table> -->
<table>
  <tr>
  <h2>human_ensembl_msigdb_c2<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c2.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c2.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c2.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c2.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c2.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c2.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_c3<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c3.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c3.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c3.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c3.R.png width='500'>
    6h_delta_exp
<!--     <td><img src=72h_delta_exp/human_ensembl_msigdb_c3.R.png width='500'>
    72h_delta_exp -->
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c3.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_c4<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c4.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c4.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c4.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c4.R.png width='500'>
    6h_delta_exp
<!--     <td><img src=72h_delta_exp/human_ensembl_msigdb_c4.R.png width='500'>
    72h_delta_exp -->
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c4.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_c5<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c5.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c5.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c5.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c5.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c5.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c5.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_c6<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c6.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c6.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c6.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c6.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c6.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c6.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_c7<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c7.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c7.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c7.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c7.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c7.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c7.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_full<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_full.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_full.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_full.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_full.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_full.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_full.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_msigdb_h<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_h.L.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_h.L.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_h.L.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_h.R.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_msigdb_h.R.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_msigdb_h.R.png width='500'>
    120h_delta_exp
  <tr>
<table>

<table>
  <tr>
  <h2>human_ensembl_RBPs_all_gene_ids<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_all_gene_ids.all.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_RBPs_all_gene_ids.all.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_RBPs_all_gene_ids.all.png width='500'>
    120h_delta_exp
  <tr>
<table>

<!-- <table>
  <tr>
  <h2>human_ensembl_RBPs_all_gene_names<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_all_gene_names.all.png width='500'>
    <td><img src=72h_delta_exp/human_ensembl_RBPs_all_gene_names.all.png width='500'>
    <td><img src=120h_delta_exp/human_ensembl_RBPs_all_gene_names.all.png width='500'>
  <tr>
<table> -->
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids.all.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids.all.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids.all.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_3UTR<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_5UTR<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png width='500'>
    120h_delta_exp
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_coding_exons<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png width='500'>
    6h_delta_exp
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png width='500'>
    72h_delta_exp
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png width='500'>
    120h_delta_exp
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_introns<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png width='500'>
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png width='500'>
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png width='500'>
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_DeepBind<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_DeepBind.all.png width='500'>
    <td><img src=72h_delta_exp/human_ensembl_RBPs_DeepBind.all.png width='500'>
    <td><img src=120h_delta_exp/human_ensembl_RBPs_DeepBind.all.png width='500'>
  <tr>
<table>
 -->
