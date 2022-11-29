inDir=$1
exp=$2
geneset=$3

cd $inDir

outDir=${exp/.txt/}; 
outDir=${outDir}_onePAGE_${geneset} 

export TEISERDIR='/data_gilbert/home/aarab/Workflows/TEISERv1.1'
echo $exp $geneset 
base=`basename $exp` 
base=${base/.txt/} 


perl ${TEISERDIR}/run_mi_gene_list.pl \
    --expfile=$exp \
    --genefile=${geneset}.txt \
    --exptype=continuous \
    --ebins=11 \
    --draw_min=-3 \
    --draw_max=3 \
    --species=human \
    --doremovedups=0 \
    --doremoveextra=0

rm -r $outDir 
mv ${exp}_GENESET $outDir 
mv ${outDir}/${base}.txt.summary.pdf ${outDir}.summary.pdf

echo 'done!'