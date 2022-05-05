PDIR=$1
fastqDIR=$2
qcDIR=$3
JOBS=$4

cd $PDIR

mkdir -p ${qcDIR}/;

for fq_file in ${fastqDIR}/*.fastq.gz; do
    fq_base=`basename $fq_file`;
    sample_id=${fq_base/.fastq.gz/}
    echo '--------------' $sample_id '--------------'

    fastqc -q -t $JOBS -o ${qcDIR}/ $fq_file;
done

multiqc ${qcDIR}/ -n mutiqc-${qcDIR}
