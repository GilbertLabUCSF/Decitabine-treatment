PDIR=$1
FASTQDIR=$2
INDEX=$3
JOBS=$4

cd $PDIR
mkdir -p ./quants/

for f in $FASTQDIR/*fastq.gz; do 
	samp=`basename ${f}`; 
	samp=${samp/.fastq.gz/}; 
	echo "Processing sample ${samp}"; 
	salmon quant -i $INDEX \
 	-l A -r $f -p $JOBS --validateMappings -o ./quants/$samp; 
done
