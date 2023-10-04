#!/bin/bash -e
###
 # @Authour: utao.cao
### 

source ~/software/miniconda3/bin/activate
conda activate ~/conda_env/atac_wkflo/

export SOURCE_PATH="./wd_cyt/sheep/ATAC"
DEEPTOOLS_DIR=${SOURCE_PATH}/deeptools
CORE_NUM=32
CORE_FUL=64
export PATH="./software/ucscTools:$PATH"
export SAMTOOLS=./conda_env/wgbs_wkflo/bin/samtools
export SAMBAMBA=./software/sambamba-0.8.0
export GENOME_PATH="./reference/sheep"

ReadsFilter=${SOURCE_PATH}/scripts/reads_filter.pl
refGENOME="${GENOME_PATH}/ncbi_sheep_2015/GCA_000298735.2_Oar_v4.0_genomic.fa.gz"
metaInfoGENOME=${SOURCE_PATH}/NCBI2015_SheepGenomeMetaInfo.txt

cd $SOURCE_PATH

for id in $(ls ${SOURCE_PATH}/align/*_sortDedup.bam)  ;do 
{
 sampleName=${id##*/} ; echo $sampleName
 bedtools bamtobed -i $id | sort -k4,4V >  ${SOURCE_PATH}/align/bedFilter/${sampleName%%_*}.sort.bed   
 perl $ReadsFilter -i ${SOURCE_PATH}/align/bedFilter/${sampleName%%_*}.sort.bed -o ${SOURCE_PATH}/align/bedFilter/${sampleName%%_*}.filter.bed -len 150
 bedGraphToBigWig  ${SOURCE_PATH}/align/bedFilter/${sampleName%%_*}.filter.bed $metaInfoGENOME ${DEEPTOOLS_DIR}/${sampleName%%_*}_notcall.bw
} &
 done



computeMatrix reference-point \
 -S `ls ${DEEPTOOLS_DIR}/*_notcall.bw | xargs echo` \
 -R $REF_UCSC \
 --referencePoint TSS \
 -a 3000 -b 3000 \
 -out ${DEEPTOOLS_DIR}/matrixATAC_compare4stages_TSS.tab.gz  &
 
computeMatrix scale-regions \
  -R $REF_UCSC \
  -S ${DEEPTOOLS_DIR}/*_notcall.bw  \
  -b 3000 -a 3000 \
  --regionBodyLength 5000 \
  --skipZeros -o ${DEEPTOOLS_DIR}/matrixATAC_compare4stages_scaledRegion.tab.gz  \
  --outFileNameMatrix ${DEEPTOOLS_DIR}/matrixATAC_compare4stages_scaledRegion.tab \
  --outFileSortedRegions ${DEEPTOOLS_DIR}/regionsGenes_compare4stages_scaledRegion.bed  &

wait
 for matrix_file in $(/bin/ls $DEEPTOOLS_DIR/matrix*tab.gz ); do
 profile_file=profile${matrix_file#*matrix}
 if [ ! -f ${profile_file%tab.gz*}pdf ];then
 echo ${profile_file%tab.gz*}pdf
 plotProfile -m $matrix_file  --perGroup -out ${profile_file%tab.gz*}pdf \
	--startLabel '' --endLabel '' --plotFileFormat  pdf \
	--samplesLabel D3 M12 M3 M6  &
 sleep 3s
 fi
  done

 for matrix_file in $(/bin/ls $DEEPTOOLS_DIR/matrix*tab.gz ); do
 profile_file=profile${matrix_file#*matrix}
 if [ ! -f ${profile_file%tab.gz*}png ];then
 echo ${profile_file%tab.gz*}png
 plotProfile -m $matrix_file  --perGroup -out ${profile_file%tab.gz*}png \
	--startLabel '' --endLabel ''   --plotFileFormat  png \
	--samplesLabel D3 M12 M3 M6 &
 sleep 3s
 fi
  done
wait
