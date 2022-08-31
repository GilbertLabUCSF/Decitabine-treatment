# https://github.com/Kingsford-Group/scallop
mkdir -p scallop/tinat/compare

# export PATH='/sadra/goodarzilab/abe/Workflows/rnaseqtools-1.0.3/bin/':$PATH

conda activate alignment
cd scallop-genome/

### 0. download genomes
# # mkdir -p hg19
# # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.basic.annotation.gtf.gz
gunzip gencode.v34lift37.basic.annotation.gtf.gz
# # mv hg19.knownGene.gtf.gz hg19/hg19.knownGene.gtf.gz
# # gunzip hg19/hg19.knownGene.gtf.gz

### 1. merge
gtfmerge union gtf-list.txt DAC_GRCh37_merged.gtf -t 2

# ### 2. compare 
# gffcompare -o hg19_tinat/gffall -r /rumi/shams/genomes/hg19/hg19_genes.gtf tinat/hg19.tinat.gtf 

# ### 3. subset
# gtfcuff puniq \
#     tinat/gffall.hg19.tinat.gtf.tmap \
#     tinat/hg19.tinat.gtf \
#     /rumi/shams/genomes/hg19/hg19_genes.gtf \
#     tinat/unique.gtf

### 4. gtf2fasta
# gffread tinat/unique.gtf -g /rumi/shams/genomes/hg19/hg19.fa -w tinat/unique.fa
gffread DAC_GRCh37_merged.gtf -g /rumi/shams/genomes/hg19/hg19.fa -w DAC_GRCh37_merged.fa

# # 5.concatenate with hg19 fasta
# gffread hg19/hg19.knownGene.gtf -g /rumi/shams/genomes/hg19/hg19.fa -w hg19/hg19.knownGene.fa
# cat hg19/hg19.knownGene.fa tinat/unique.fa > tinat/hg19.tinat.fa

### aligner index 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip GRCh37.primary_assembly.genome.fa.gz

mkdir STAR.index

STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir STAR.index \
    --genomeFastaFiles GRCh37.primary_assembly.genome.fa \
    --sjdbGTFfile DAC_GRCh37_merged.gtf \
    --sjdbOverhang 99

rm GRCh37.primary_assembly.genome.fa.gz
