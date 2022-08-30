mkdir -p scallop/tinat/compare

PATH=$PATH:/rumi/shams/abe/Workflows/rnaseqtools-1.0.3/bin/

conda activate alignment

# 1. merge
gtfmerge union scallop/tinat/tinat-gtf-list.txt scallop/tinat/hg19.tinat.gtf -t 18
# 2. compare 
gffcompare -o scallop/tinat/gffall -r /rumi/shams/genomes/hg19/hg19_genes.gtf scallop/tinat/hg19.tinat.gtf 
# 3. subset
gtfcuff puniq scallop/tinat/gffall.hg19.tinat.gtf.tmap scallop/tinat/hg19.tinat.gtf /rumi/shams/genomes/hg19/hg19_genes.gtf scallop/tinat/unique.gtf
# 4. gtf2fasta
gffread scallop/tinat/unique.gtf -g /rumi/shams/genomes/hg19/hg19.fa -w scallop/tinat/unique.fa
# 5.concatenate with hg19 fasta
# mkdir -p scallop/hg19
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz
# mv hg19.knownGene.gtf.gz scallop/hg19/hg19.knownGene.gtf.gz
# gunzip scallop/hg19/hg19.knownGene.gtf.gz
gffread scallop/hg19/hg19.knownGene.gtf -g /rumi/shams/genomes/hg19/hg19.fa -w scallop/hg19/hg19.knownGene.fa
cat scallop/hg19/hg19.knownGene.fa scallop/tinat/unique.fa > scallop/tinat/hg19.tinat.fa
