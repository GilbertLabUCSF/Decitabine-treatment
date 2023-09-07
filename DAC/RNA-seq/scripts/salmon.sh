PDIR=$1
FASTQDIR=$2
INDEX=$3
JOBS=$4

quantsDir='exp/quants/'

cd $PDIR
mkdir -p $quantsDir

for f in $FASTQDIR/*fastq.gz; do 
	samp=`basename ${f}`; 
	samp=${samp/.fastq.gz/}; 
	echo "Processing sample ${samp}"; 
	salmon quant -i $INDEX \
 	-l A -r $f -p $JOBS --validateMappings -o $quantsDir/$samp; 
done
