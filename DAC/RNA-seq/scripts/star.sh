index=$1
fastqDIR=$2
bamDIR=$3
JOBS=$4

mkdir -p ${bamDIR}
mkdir -p ${bamDIR}_star_qc

STAR --genomeLoad LoadAndExit --genomeDir $index

for fq in ${fastqDIR}/*fastq.gz; do
    fq=`basename $fq`
    out=${fq/.fastq.gz/}
    echo '------------' $out '-----------'
    STAR \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --runThreadN $JOBS \
    --genomeDir $index \
    --readFilesIn ${fastqDIR}/$fq \
    --outFileNamePrefix ${bamDIR}/$out
    
    mv -v ${bamDIR}/${out}Aligned.sortedByCoord.out.bam ${bamDIR}/${out}.bam
    mv -v ${bamDIR}/${out}Log.final.out ${bamDIR}_star_qc/
    rm -v ${bamDIR}/${out}*out*
    rm -r ${bamDIR}/${out}_STARtmp/
    
done

STAR --genomeLoad Remove --genomeDir $index

rm -r _STARtmp/ Log.out Log.progress.out Aligned.out.sam
