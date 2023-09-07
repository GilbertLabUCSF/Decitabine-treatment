index=$1
bamDIR=$2
JOBS=$3

mkdir -p ${bamDIR}
mkdir -p ${bamDIR}_star_qc

STAR --genomeLoad LoadAndExit --genomeDir $index

for fq in fastq/*R1*; do
    fq=`basename $fq`
    out=${fq/_R1*/}
    echo '------------' $out '-----------'
    STAR \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --runThreadN $JOBS \
    --genomeDir $index \
    --readFilesIn fastq/$fq \
    --outFileNamePrefix ${bamDIR}/$out
    
    mv -v ${bamDIR}/${out}Aligned.sortedByCoord.out.bam ${bamDIR}/${out}.bam
    mv -v ${bamDIR}/${out}Log.final.out ${bamDIR}_star_qc/
    rm -v ${bamDIR}/${out}*out*
    rm -rv ${bamDIR}/${out}_STARtmp/
    
done

STAR --genomeLoad Remove --genomeDir $index

rm -r _STARtmp/ Log.out Log.progress.out Aligned.out.sam
