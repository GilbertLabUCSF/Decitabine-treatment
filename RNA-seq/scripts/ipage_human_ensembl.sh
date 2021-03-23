exp_file=$1;
outdir=${exp_file/.txt/_dependent}

ipage_ann='/flash/bin/iPAGEv1.0/PAGE_DATA/ANNOTATIONS/'

mkdir -p $outdir

for f in `ls -d ${ipage_ann}/human_ensembl*`; do

    base=`basename "$f"`;
    echo '________________' $base '________________';

    if [ -d "${outdir}/${base}/" ];
    then
        echo 'This result exist!';
    ## TO-DO; Force run to remove stuff from previous run.
    else
        # Run iPAGE
        perl $PAGEDIR/page.pl --expfile=$exp_file --independence=0 \
        --species=$base --exptype=continuous --ebins=11 --nodups=1; 
        wait
        mv -v ${exp_file}_PAGE/ ${outdir}/${base}/;

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
            if [ "$(wc -l < ${outdir}/${base}/pvmatrix.${side}.txt)" -eq "$(echo '1')" ];
            then
                echo '_____...________' $base $side '-> No signal!';
            else
                echo '_____...________' $base $side '...';
                mkdir -p ${exp_file}_PAGE;
                perl /flash/bin/iPAGEv1.0/SCRIPTS/mi_go_draw_matrix.pl  \
                --pvaluematrixfile=${outdir}/${base}/pvmatrix.${side}.txt \
                --expfile=${exp_file} \
                --order=1 --draw_sample_heatmap=false \
                --min=-3 --max=3 \
                --cluster=5 --quantized=0;
                wait
                for sum in `ls ${exp_file}_PAGE`; do
                    sumS=${sum/summary/$side};
                    mv -v ${exp_file}_PAGE/${sum} ${outdir}/${base}/${sumS};
                done
                rm -vr ${exp_file}_PAGE;
            fi;

        done

    fi


    # pdf2png!
    for file in ${outdir}/${base}/${exp_file}*.pdf; do 
        file=${file/.pdf/}
        bash /rumi/shams/abe/Workflows/my_scripts/pdf2png.sh ${file}.pdf
        new=${file/\/$exp_file/}
        cp -v ${file}.pdf ${new}.pdf
        cp -v ${file}.png ${new}.png
    done
done

for pv in `ls ${outdir}/*/pvmatrix.*.txt`; do
    if  [ "$(wc -l < $pv)" -eq "$(echo '1')" ]; then
        echo `dirname $pv` $b NO SIGNAL!
        rm -v $pv
    elif [ "$(wc -l < $pv)" -eq "$(echo '0')" ]; then
        echo `dirname $pv` $b NO DATA!
        rm -v $pv
    else
        echo `dirname $pv` $b `cat $pv | wc -l`
    fi
done

echo `date`
echo "all done!"
