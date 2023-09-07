PDIR=$1
JOBS=$2
INDEX=$3

cd $PDIR
mkdir -p bam
mkdir -p logs/star_aligner


STAR --genomeLoad LoadAndExit --genomeDir $INDEX
for fq in fastq/*R1*fastq.gz; do
    fq1=`basename $fq`
    fq2=${fq1/_R1_/_R2_}
    out=${fq1/_R1_001.fastq.gz/}
    echo '--------------' $out '--------------'
    STAR \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --runThreadN $JOBS \
    --genomeDir $INDEX \
    --readFilesIn fastq/$fq1 fastq/$fq2 \
    --outFileNamePrefix bam/$out \
    --outReadsUnmapped Fastx;

    mv -v bam/${out}Aligned.sortedByCoord.out.bam bam/${out}.bam
    mv -v bam/${out}Log.final.out logs/star_aligner/
    rm -v bam/${out}*out*
done

STAR --genomeLoad Remove --genomeDir $INDEX

rm -r _STARtmp/ Log.out Log.progress.out Aligned.out.sam
