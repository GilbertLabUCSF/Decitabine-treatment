work_dir=$1
cd $work_dir

echo "step 1: prepare inputs for FIRE"
python $TEISERDIR/prep_fasta_for_fire_run.py peak.fa
perl $TEISERDIR/prep_seqs_for_teiser_run.pl peak.fa peaks
echo "--- DONE! ---"

echo "step 2: run FIRE for known m6A motifs (non-discovery mode)"
perl $FIREDIR/fire.pl --expfile=peak_fire.txt --exptype=discrete --fastafile_rna=peak_fire.fa \
--nodups=1 --dodna=0 --dodnarna=0 --species=human --doskipdiscovery=1 \
--motiffile_rna=$MOTIF --oribiasonly=0 > non-discovery_FIRE.log
rm -rv non-discovery_FIRE 
mv -v peak_fire.txt_FIRE/ non-discovery_FIRE 
echo "--- DONE! ---"

# echo "step 3: run FIRE discovery mode"
# perl $FIREDIR/fire.pl --expfile=peak_fire.txt --exptype=discrete --fastafile_rna=peak_fire.fa \
# --nodups=1 --dodna=0 --dodnarna=0 --species=human --oribiasonly=0 > discovery_FIRE.log 
# rm -rv discovery_FIRE 
# mv -v peak_fire.txt_FIRE/ discovery_FIRE
# echo "--- DONE! ---"
