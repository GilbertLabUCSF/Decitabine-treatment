FASTQDIR=$1
OUTDIR=$2

INDEX=/rumi/shams/abe/Gilbertlab/Decitabine_treatment/RNA-seq/hl60-exp/scallop/tinat/salmon.index

mkdir -p $OUTDIR

for f in $FASTQDIR/*fastq.gz; do 
	samp=`basename ${f}`; 
	samp=${samp/.fastq.gz/}; 
	echo "Processing sample ${samp}"; 
	salmon quant -i $INDEX \
 	-l A -r $f -p 18 -o ${OUTDIR}/${samp}; 
done
