#!/bin/sh
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -J bsPMD.cov
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log

covlist=$(ls ../*bismark.cov | tr '\n' " ")
path=/public4/home/sc30941/wd_cyt/ay_sheep/pj013_WGBS/align_bismark/unsort_Bismarkfile/methylation_extraction/unionMethyl
#for covfile in $covlist;do
#  cov=${covfile#.*/}
#	awk -v OFS="\t" -v c4=${cov%%_*} '{if(NR!=1 && $1 !~ "A|J" ) {$1="chr"$1; print($1,$2,$3,$4);}else if(NR == 1) {$1="chr"$1; print "chr\tstart\tstop",c4"\n"$1,$2,$3,$4 }}' \
#          ../$cov >  ${cov%%_*}_bsmkcov_simple.bed  &
# done	 

# after obtain mmerged file
cd $path
for covfile in $covlist;do
    cov=${covfile#.*/}
    # c=0 help to calculate $3=N $4=X
	cat ../$cov | awk  -v OFS='\t' -f /public4/home/sc30941/wd_cyt/ay_sheep/pj013_ATAC/script_revertChromosome.awk   | cut -f 1,2,5,6 | \
					awk -v OFS="\t" -v c=0  'NR==FNR{c=$3;$3=$3+$4;$4=c;a[$1,$2]=$0;next}NR>FNR{if($1"\034"$2 in a)print a[$1,$2]}'  - $path/bsmkcov_AllSample_merged_position.txt >  $path/${cov%%_*}_bs_PMD.cov  && free -h &
    
done
wait


