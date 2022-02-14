## $1: input file contain table with first two columns as gene name/id and numeric value
pattern='msigdb*'
# pattern='human_*_gs*'

expfile=`basename $1`
outdir=${1/.txt/};

mkdir -p $outdir
cd $outdir; cd ../

for f in `ls -d $PAGEDIR/PAGE_DATA/ANNOTATIONS/${pattern}`; do

    base=`basename "$f"`;
    echo '________________' $base '________________';

    if [ -d "${outdir}/${base}/" ]; 
    then
        echo 'This result exist!';
    else
        # Run iPAGE 
        perl $PAGEDIR/page.pl --expfile=$expfile --species=$base --exptype=continuous --ebins=11 --nodups=1;
        # --independence=0; option for comparing results between multiple smaples.
        wait

        mv -v ${expfile}_PAGE/ ${outdir}/${base}/;

        # remove the result folder if it was empty
        counter="$(wc -l < ${outdir}/${base}/pvmatrix.txt)"
        if [ $counter -le "$(echo '1')" ] || [ -z $counter ]
        then
            rm -r ${outdir}/${base}/
        else
            # keep complete results 
            pv=${outdir}/${base}/pvmatrix.txt;
            pv0=${pv/.txt/.all.txt}; pvL=${pv/.txt/.L.txt}; pvR=${pv/.txt/.R.txt};

            # subset pathways with p-value > 2 in the first (Left) or last (Right) cluster
            mv -v $pv $pv0;
            cat $pv0 | awk -F'\t' 'BEGIN{FS=OFS="\t"}; $2>2 {print $0}' > $pvL;
            cat $pv0 | awk -F'\t' 'BEGIN{FS=OFS="\t"}; $12>2 {print $0}' > $pvR;

            for sum in `ls ${outdir}/${base}/*.summary*`; do 
                sum0=${sum/summary/all};
                mv -v $sum $sum0;
            done
            wait

            # draw sided heatmaps     
            declare -a Sides=('L' 'R');

            for side in "${Sides[@]}"; do

                # include if pv matrix contain any values 
                counter="$(wc -l < ${outdir}/${base}/pvmatrix.${side}.txt)"
                if [ $counter -le "$(echo '1')" ] || [ -z $counter ]
                then 
                    echo '_____...________' $base $side '-> No signal!';    
                else 
                    echo '_____...________' $base $side '...';
                    mkdir -p ${expfile}_PAGE;
                    perl $PAGEDIR/SCRIPTS/mi_go_draw_matrix.pl  \
                    --pvaluematrixfile=${outdir}/${base}/pvmatrix.${side}.txt \
                    --expfile=$expfile \
                    --order=1 --draw_sample_heatmap=false \
                    --min=-3 --max=3 \
                    --cluster=5 --quantized=0;
                    wait
                    for sum in `ls ${expfile}_PAGE/*.summary*.pdf`; do 
                        mv -v ${sum} ${outdir}/${base}.${side}.pdf
                    done
                    rm -vr ${expfile}_PAGE;
                fi;

            done
            
            mv -v $pv0 $pv
            
        fi;
        
    fi

done
