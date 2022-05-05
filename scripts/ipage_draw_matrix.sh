function draw_matrix()
{
    export PAGEDIR='/data_gilbert/home/aarab/iPAGE'
    
    mkdir -p ${1}_PAGE;

    local expfile=$1;
    local pvaluematrixfile=$2;
    
    perl ${PAGEDIR}/SCRIPTS/mi_go_draw_matrix.pl  \
            --pvaluematrixfile=${pvaluematrixfile} \
            --expfile=${expfile} \
            --order=1 --draw_sample_heatmap=false \
            --min=-3 --max=3 \
            --cluster=5 --quantized=0;
    wait
    pic=${expfile/.txt/.pdf}
    for sum in `ls ${expfile}_PAGE/*.summary*.pdf`; do 
        mv -v ${sum} ${pic}
    done
    rm -vr ${expfile}_PAGE;
}

draw_matrix $1 $2