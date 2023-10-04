#!/bin/bash
###
 # @Authour: utao.cao
### 
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -J pivot
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log

# example:
#zcat pmdFullAnnoStat_100bp_4category.txt.gz | grep --color=auto -e "^chr21" | cut -f 1-15 > wgbs_meanMethy100bp_chr21.txt 

cd ./wd_cyt/.sheep/pj013_WGBS/align_bismark/unsort_Bismarkfile/methylation_extraction/unionMethyl/pmd_ReanalysisResult
python  ./parallel/dissectDateFrame4parallel.py wgbs_meanMethy100bp_chr22.txt  __TMP  256
cd __TMP
parallel -k python ../parallel/pivotSummary_parallelVersion.py {}  output ::: `ls | grep tmp`
python ../parallel/pivotTableMerge.py ./output/  DissectedDF
rm tmpDiseectedDF_*
wait
