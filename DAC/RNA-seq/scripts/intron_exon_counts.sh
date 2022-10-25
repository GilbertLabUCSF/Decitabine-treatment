#loop featureCounts intron and exon counts

PDIR=$1
bamDIR=$2
countDIR=$3
logDIR=$4
JOBS=$5

cd $PDIR

GTF_index='/data_gilbert/home/aarab/genomes/hg38/hg38_ensemble_'

mkdir -p $countDIR
mkdir -p $logDIR

for f in ${bamDIR}/*.bam; do
    base=`basename $f`
    sample=${base/.bam/};
    echo -e '----------------------- ' $sample  ' -----------------------'
    echo `date` 
    featureCounts -M -T $JOBS -t intron -g gene_id -a ${GTF_index}introns.gtf -o ${countDIR}/${sample}_introns.txt ${f} &> ${logDIR}/${sample}_introns.log
    featureCounts -M -T $JOBS -t exon -g gene_id -a ${GTF_index}consExons.gtf -o ${countDIR}/${sample}_exons.txt ${f} &> ${logDIR}/${sample}_exons.log
    echo `date` Done!
done
