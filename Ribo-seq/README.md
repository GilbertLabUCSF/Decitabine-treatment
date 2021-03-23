Tools: 
- https://github.com/goodarzilab/Ribolog
- https://github.com/goodarzilab/PAGE

```bash
cat fit1_deci_fdr_qval_g_0.csv | awk -F, '{print $2"\t"$4}' > lnTE_T_vs_U.txt
```

# QC plots

<table>
    <tr>
    <td><img src="plots/QC_and_volcano_plots-0.png" style="width:500px">
    <td><img src="plots/QC_and_volcano_plots-1.png" style="width:500px">
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/QC_and_volcano_plots-2.png" style="width:500px">
    <td><img src="plots/QC_and_volcano_plots-3.png" style="width:500px">
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/QC_and_volcano_plots-4.png" style="width:500px">
    <td><img src="plots/QC_and_volcano_plots-5.png" style="width:500px">
    <tr>
<table>

<h2>DMSO, rep 1<h2>
<table>
    <tr>
    <td><img src="plots/deci_RPF_read_end_heatmaps-0.png" style="width:600px">
    RPF_read_end_heatmaps
    <td><img src="plots/deci_RPF_Read_length_distributions-0.png" style="width:300px">
    RPF_Read_length_distributions
    <tr>
<table>

<table>
    <tr>
    <td><img src="plots/Periodicity_by_length_region-0.png" style="width:500px">
    RPF_read_end_heatmaps
    <td><img src="plots/Periodicity_by_region-0.png" style="width:400px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/RPF_ribosome_occupancy_profiles_from_annotation-0.png" style="width:600px">
    RPF_ribosome_occupancy_profiles_from_annotation
    <tr>
<table>

<h2>DMSO, rep 2<h2>
<table>
    <tr>
    <td><img src="plots/deci_RPF_read_end_heatmaps-1.png" style="width:600px">
    RPF_read_end_heatmaps
    <td><img src="plots/deci_RPF_Read_length_distributions-1.png" style="width:300px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/Periodicity_by_length_region-1.png" style="width:500px">
    RPF_read_end_heatmaps
    <td><img src="plots/Periodicity_by_region-1.png" style="width:400px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/RPF_ribosome_occupancy_profiles_from_annotation-1.png" style="width:600px">
    RPF_ribosome_occupancy_profiles_from_annotation
    <tr>
<table>

<h2>Drug, rep 1<h2>
<table>
    <tr>
    <td><img src="plots/deci_RPF_read_end_heatmaps-2.png" style="width:600px">
    RPF_read_end_heatmaps
    <td><img src="plots/deci_RPF_Read_length_distributions-2.png" style="width:300px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/Periodicity_by_length_region-2.png" style="width:500px">
    RPF_read_end_heatmaps
    <td><img src="plots/Periodicity_by_region-2.png" style="width:400px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/RPF_ribosome_occupancy_profiles_from_annotation-2.png" style="width:600px">
    RPF_ribosome_occupancy_profiles_from_annotation
    <tr>
<table>

<h2>Drug, rep 2<h2>
<table>
    <tr>
    <td><img src="plots/deci_RPF_read_end_heatmaps-3.png" style="width:600px">
    RPF_read_end_heatmaps
    <td><img src="plots/deci_RPF_Read_length_distributions-3.png" style="width:300px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/Periodicity_by_length_region-3.png" style="width:500px">
    RPF_read_end_heatmaps
    <td><img src="plots/Periodicity_by_region-3.png" style="width:400px">
    RPF_Read_length_distributions
    <tr>
<table>
<table>
    <tr>
    <td><img src="plots/RPF_ribosome_occupancy_profiles_from_annotation-3.png" style="width:600px">
    RPF_ribosome_occupancy_profiles_from_annotation
    <tr>
<table>

# Enrichment analysis 
<table>
  <tr>
  <h2>human_ensembl<h2>
  <img src=lnTE_T_vs_U/human_ensembl.L.png style="width:600px">
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_encode_tf<h2>
  <img src=lnTE_T_vs_U/human_ensembl_encode_tf.all.png style="width:600px">
  <tr>
<table> -->
<!-- <table>
  <tr>
  <h2>human_ensembl_msigdb_c1<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c1.all.png style="width:600px">
  <tr>
<table> -->
<table>
  <tr>
  <h2>human_ensembl_msigdb_c2<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c2.R.png style="width:600px">
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c2.L.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c3<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c3.L.png style="width:600px">
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_msigdb_c4<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c4.all.png style="width:600px">
  <tr>
<table> -->
<table>
  <tr>
  <h2>human_ensembl_msigdb_c5<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c5.all.png style="width:600px">
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c5.L.png style="width:600px">
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_msigdb_c6<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c6.all.png style="width:600px">
  <tr>
<table> -->
<table>
  <tr>
  <h2>human_ensembl_msigdb_c7<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c7.all.png style="width:600px">
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_c7.L.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_full<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_full.R.png style="width:600px">
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_full.L.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_h<h2>
  <img src=lnTE_T_vs_U/human_ensembl_msigdb_h.all.png style="width:600px">
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_RBPs_all_gene_ids<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_all_gene_ids.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_all_gene_names<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_all_gene_names.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_coding_gene_ids.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_3UTR<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png style="width:600px">
  <tr>
<table> -->
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_5UTR<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png style="width:600px">
  <tr>
<table>
<!-- <table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_coding_exons<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_introns<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_DeepBind<h2>
  <img src=lnTE_T_vs_U/human_ensembl_RBPs_DeepBind.all.png style="width:600px">
  <tr>
<table>
 -->
