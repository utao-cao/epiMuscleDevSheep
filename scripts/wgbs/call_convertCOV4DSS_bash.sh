#!/bin/sh
###
 # @Authour: utao.cao
### 

pyS='./wd_cyt/sheep/pj013_WGBS/script/convertCOV4DSS.py'
covD='./wd_cyt/sheep/pj013_WGBS/align_bismark/unsort_Bismarkfile/methylation_extraction'
metaInfoGENOME="./wd_cyt/sheep/pj013_ATAC/script_revertChromosome.awk"
depthThreshold=5

# ɸѡcovdepth, select autosome, control %.2f
for covF in $(ls $covD/*bismark.cov);do
   echo "deal $covF"
   basename=${covF##*/}
   # 1. Convert 4 DSS
    #python3 $pyS $covF 
    # 2. prepare for annotation statistical
    cat $covF | awk  -v OFS='\t'  -f  $metaInfoGENOME | \
      awk  -v OFS='\t' -v depthT="$depthThreshold" \
            '{if( $5+$6 > depthT ) printf ("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$4)}' >  ${pyS%/*}/DSS/${basename%%_*}_bsmk_ovis.cov
done

# insert header may cause error
# cat ${covF}  | awk -v OFS="\t"  -f  $metaInfoGENOME - | sed '1i\chr\tpos\tN\tX'   >  ${pyS%/*}/DSS/${basename%%_*}_4DSS_ovis.cov
