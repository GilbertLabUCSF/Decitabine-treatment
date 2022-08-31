index=$1
bamDIR=$2

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
    --runThreadN 18 \
    --genomeDir $index \
    --readFilesIn fastq/$fq \
    --outFileNamePrefix ${bamDIR}/$out
done

STAR --genomeLoad Remove --genomeDir $index

rm -vr _STARtmp/ Aligned.out.sam Log.out Log.progress.out
