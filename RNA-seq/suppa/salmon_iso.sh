PDIR=$1
quants=$2
JOBS=$3

FASTQDIR='fastq'
INDEX='/data_gilbert/home/aarab/genomes/hg38/gencode.v34/salmon_index/'

cd $PDIR
mkdir -p $quants

for fq in $FASTQDIR/*.fastq.gz; do
	samp=`basename ${fq}`;
    samp=${samp/.fastq.gz/};
    # fq2=${fq1/_R1/_R2};
	echo "Processing sample ${samp}";
	cmd="salmon quant -i $INDEX -l ISF --gcBias -r $fq -p $JOBS -o $quants/$samp"
	echo $cmd
	$cmd &> $quants/${samp}.log
	echo DONE@ `date`;
done
