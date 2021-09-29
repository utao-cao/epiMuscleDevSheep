#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -J pivot
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
cd /public4/home/sc30941/wd_cyt/ay_sheep/pj013_WGBS/align_bismark/unsort_Bismarkfile/methylation_extraction/unionMethyl/pmd_ReanalysisResult
python  ./parallel/diseectDateFrame4parallel.py pmdFullAnnoStat_100bp_4category.txt  __TMP  256
cd __TMP
parallel -k python ../parallel/pivotSummary_parallelVersion.py {}  output ::: `ls | grep tmp`
python ../parallel/pivotTableMerge.py ./output/  DissectedDF
rm tmpDiseectedDF_*
wait