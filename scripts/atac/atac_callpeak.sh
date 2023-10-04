#!/bin/sh
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -J atacalPk
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
## while can not parrallel, so replaced with parrallel
#links:https://stackoverflow.com/questions/30760449/gnu-parallel-redirect-output-to-a-file-with-a-specific-name   answer02
# find ${SOURCE_PATH}/peaks/ -type f  -name "*_sortbdg.bed" | parallel -j 6 "awk -f  ${SOURCE_PATH}/script.awk {} " >  '>' {/.}.txt

set -e
export SOURCE_PATH="./wd_cyt/sheep/ATAC"
CORE_NUM=32
CORE_FUL=64
export PATH="./software/ucscTools:$PATH"
export SAMTOOLS=./conda_env/wgbs_wkflo/bin/samtools
export SAMBAMBA=./software/sambamba-0.8.0
export GENOME_PATH="./reference/sheep"

###
# change
###
## ���Ŀ¼��Ҫ����
ReadsFilter=${SOURCE_PATH}/scripts/reads_filter.pl
BT2_INDEX=${SOURCE_PATH}/reference/Sheep_bowtie
refGENOME=${GENOME_PATH}/Ovis_aries.Oar_v3.1.dna.toplevel.fa.gz
metaInfoGenme=${SOURCE_PATH}/SheepGenome_4.0.2015_MetaInfo.txt
REF_UCSC=${SOURCE_PATH}/reference/ucsc.Husheep.bed

metaInfoGENOME=${SOURCE_PATH}/HuSheepGenomeMetaInfo.txt
cd ${SOURCE_PATH}
sh ${SOURCE_PATH}/clipBED_autoGenerateAWKScript.sh $metaInfoGENOME #SheepGenome_4.0.2015_MetaInfo.txt

GSIZE=$(bowtie2-inspect $BT2_INDEX |perl -ne 'BEGIN{$n=0} next if(/^>/);s/[^ATGC]//gi; $n+=length($_); END{print int($n*0.85);}')
echo "Effective Genome Size" $GSIZE

DIRLISTs='align/bedFilter clean reference tss/Body deeptools_result'
for DIR_LIST in $DIRLISTs ;do
    if [ ! -d $DIR_LIST ];then
      mkdir -p $DIR_LIST
    fi
done

##########################
# 02_index
##########################
if [ ! -f ${BT2_INDEX}* ];then
 bowtie2-build -f $refGENOME --threads $CORE_FUL  ${BT2_INDEX}
fi

##########################
# 02_trim_galore.sh
##########################


ls raw/*_1.fq.gz > raw.R1
ls raw/*_2.fq.gz > raw.R2
ls raw/ | sed 's/_.*//g' | uniq > raw.0
paste raw.0 raw.R1 raw.R2 > config.raw

fastqc -t $CORE_NUM  ${SOURCE_PATH}/raw/*gz -o ${SOURCE_PATH}/raw  &

cat config.raw  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}

trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o  ${SOURCE_PATH}/clean  $fq1   $fq2  &

done

wait
ps -ef |grep trim

fastqc -t $CORE_NUM  ${SOURCE_PATH}/clean/*gz -o ${SOURCE_PATH}/clean/ &

wait

##########################
#03 align
##########################
cd ${SOURCE_PATH}/align

ls ${SOURCE_PATH}/clean/*_1.fq.gz > 1
ls ${SOURCE_PATH}/clean/*_2.fq.gz > 2
##@annotate: obtain sampleID
#@Attention: need to be modified if necessary!
ls ${SOURCE_PATH}/clean/*_2.fq.gz |cut -d"/" -f 9|cut -d"_" -f 1  > 0 ## ��������Լ��𲽷ֿ�����һ�¼��0����Ľ���Ƿ������sample����
paste 0 1 2  > config.clean ## ��mappingʹ�õ������ļ�

#source activate atac 
cat config.clean |while read id;
#@annotate: read line,and each line was constructed as a new sub-array to be index
do echo $id
   arr=($id)
   fq2=${arr[2]}
   fq1=${arr[1]}
   sample=${arr[0]}

## �ȶԹ���15����һ������
#@annotate: raw.bam -> (INDEX) raw.bed , raw.stat 
#					->  .rmdup.bam (INDEX)  ->  rmdup.stat
#					->  .last.bam(filter)   ->  .bed   .last.stat   
   echo "bowtie2"
   bowtie2  -p ${CORE_FUL}  --very-sensitive -X 2000 -x  $BT2_INDEX -1 $fq1 -2 $fq2 | ${SAMTOOLS} sort -l 9 -O bam  -@ ${CORE_FUL} -o - > ${sample}.raw.bam
   ${SAMTOOLS} index -@ ${CORE_FUL} ${sample}.raw.bam
    echo "bedtools bamtobed"
   bedtools bamtobed -i ${sample}.raw.bam  > ${sample}.raw.bed
   ${SAMTOOLS} flagstat ${sample}.raw.bam  > ${sample}.raw.stat
   
   sleep 10s
done
cd ${SOURCE_PATH}/align
cat config.clean |while read id;
do echo $id
   arr=($id)
   fq2=${arr[2]}
   fq1=${arr[1]}
   sample=${arr[0]}
echo "SAMBAMBA markdup"
   # only perform with 2 cores, and with more cores saves little time; 64 cores got stucked
   # will automaticlly index
   time $SAMBAMBA markdup  -t 8 --tmpdir='./'  -r ${sample}.raw.bam  ${sample}.rmdup.bam
   ${SAMTOOLS} index  -@ ${CORE_FUL}  ${sample}.rmdup.bam

## ref:https://www.biostars.org/p/170294/ 
## Calculate %mtDNA:
   mtReads=$(samtools idxstats  ${sample}.rmdup.bam | grep 'chrM' | cut -f 3)
   totalReads=$(samtools idxstats  ${sample}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
   echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

   ${SAMTOOLS} flagstat -@ ${CORE_FUL} ${sample}.rmdup.bam > ${sample}.rmdup.stat
   ${SAMTOOLS} view  -h  -f 2 -q 30    ${sample}.rmdup.bam   |grep -v chrM |${SAMTOOLS} sort  -O bam -l 9 -@ $CORE_NUM -o - > ${sample}.last.bam
   ${SAMTOOLS} index  -@ ${CORE_FUL} ${sample}.last.bam
   ${SAMTOOLS} flagstat  -@ ${CORE_FUL} ${sample}.last.bam > ${sample}.last.stat
   bedtools bamtobed -i ${sample}.last.bam  > ${sample}.bed
done




for id in $(ls ${SOURCE_PATH}/align/*last.bam | grep -v raw)  ;do 
{
sampleName=${id##*/} ; echo $sampleName
bedtools bamtobed -i $id | sort -k4,4V >  ${SOURCE_PATH}/align/bedFilter/${sampleName%%.*}.sort.bed   
perl $ReadsFilter -i ${SOURCE_PATH}/align/bedFilter/${sampleName%%.*}.sort.bed -o ${SOURCE_PATH}/align/bedFilter/${sampleName%%.*}.bed.filter -len 150
} &
done

wait

cd ${SOURCE_PATH}/peaks/
JOB_NAME="call peak"
echo $JOB_NAME
for id in $(ls ${SOURCE_PATH}/align/bedFilter/*bed.filter )  ;do 
{
 sampleName=${id##*/} ; echo $sampleName
 macs2 callpeak -t $id   --keep-dup=all --cutoff-analysis -g $GSIZE -B --SPMR --nomodel --shift  -75 --extsize 150  -n ${SOURCE_PATH}/peaks/${sampleName%%.*}  
 } & 
 done

wait

JOB_NAME="sort bed"
echo $JOB_NAME
ls ${SOURCE_PATH}/peaks/*.narrowPeak | grep -v _sort | while read id ;do 
 basename=${id%_*}; echo $basename "sortbdg.bed"
 sort -k1,1 -k2,2n  ${basename}_treat_pileup.bdg >  ${basename}_sortbdg.bed &
 sleep 300s
 done

wait

#create config.last_bam file
for i in $( ls ${SOURCE_PATH}/align/*bam | grep  last); do echo ${i##*/}  $i;done > ${SOURCE_PATH}/config.last_bam
#��ȡbam�ļ��ĵھ���, ����һ��config.last_bam�ļ����������ݰ���bam�ļ�������
cat ${SOURCE_PATH}/config.last_bam |while read id;
do
arr=($id)
sample=${arr[0]%%.*}
sample_name=${arr[1]##*/}
${SAMTOOLS} view $sample |awk '{print $9}'  > ${sample_name}_length.txt
done

JOB_NAME="bedGraphToBigWig"
echo $JOB_NAME
#deeptools�Ŀ��ӻ�
find ${SOURCE_PATH}/peaks/ -type f  -name "*_sortbdg.bed" | parallel --delay  6 "awk -f  ${SOURCE_PATH}/script.awk {} " '>' ${SOURCE_PATH}/peaks/{/.}.clip.bed

cd  ${SOURCE_PATH}/deeptools_result
ls ${SOURCE_PATH}/peaks/*bed | grep  _sortbdg.clip.bed |while read id;do
basename_id=${id##*/}
#using align/*last.bam
#bamCoverage -p ${CORE_FUL} --normalizeUsing RPKM -b $id -o ${SOURCE_PATH}/deeptools_result/${basename_id%%.*}.last.bw 
echo $basename_id
# while can not parrallel, so replaced with parrallel
# awk -f  ${SOURCE_PATH}/script.awk $id  > ${SOURCE_PATH}/peaks/${basename_id%%_*}.clip.bed  && \
 bedGraphToBigWig -unc ${SOURCE_PATH}/peaks/${basename_id%%.*}.clip.bed  ${metaInfoGENOME} ${SOURCE_PATH}/deeptools_result/${basename_id%%_*}.out.bw 

#rm $id
done


JOB_NAME="deeptools cal"
echo $JOB_NAME
a=(10000 3000 1500 )
basename_id="sheep"
for disTSS in ${a[@]};do
#for all samples 
#07_deeptools_TSS.sh
## both -R and -S can accept multiple files 
#-S ${SOURCE_PATH}/deeptools_result/*.last.bw  \

cd   ${SOURCE_PATH}/tss 
if [ $disTSS -eq 10000 ];then
	computeMatrix reference-point  --referencePoint TSS  -p ${CORE_FUL}  \
	-b ${disTSS} -a ${disTSS}    \
	-R $REF_UCSC  \
	-S ${SOURCE_PATH}/deeptools_result/*.out.bw  \
	--skipZeros  -o matrix${disTSS}_${basename_id%%.*}_TSS.gz  \
	--outFileSortedRegions regions${disTSS}_${basename_id%%.*}_genes.bed
 fi
##     both plotHeatmap and plotProfile will use the output from   computeMatrix
plotHeatmap -m matrix${disTSS}_${basename_id%%.*}_TSS.gz  -out ${basename_id%%.*}Tss${disTSS}Heatmap.png   &
plotHeatmap -m matrix${disTSS}_${basename_id%%.*}_TSS.gz  -out ${basename_id%%.*}Tss${disTSS}Heatmap_360.pdf --plotFileFormat pdf  --dpi 360  &
plotProfile -m matrix${disTSS}_${basename_id%%.*}_TSS.gz  -out ${basename_id%%.*}Tss${disTSS}_Profile.png   &
plotProfile -m matrix${disTSS}_${basename_id%%.*}_TSS.gz  -out ${basename_id%%.*}Tss${disTSS}_Profile_720.pdf --plotFileFormat pdf --perGroup --dpi 720 &

sleep 60s
##07_deeptools_Body.sh
#source activate atac
#-S ${SOURCE_PATH}/deeptools_result/*.last.bw  \

cd ${SOURCE_PATH}/tss/Body
if [ $disTSS -eq 10000 ];then
 computeMatrix scale-regions  -p ${CORE_FUL}  \
 -R $REF_UCSC  \
 -S ${SOURCE_PATH}/deeptools_result/*.out.bw  \
 -b ${disTSS} -a ${disTSS}  \
 --skipZeros -o matrix${disTSS}_${basename_id%%.*}_body.gz
 fi
plotHeatmap -m matrix${disTSS}_${basename_id%%.*}_body.gz  -out ${basename_id%%.*}Tss${disTSS}_body_Heatmap_360.png   &
plotHeatmap -m matrix${disTSS}_${basename_id%%.*}_body.gz  -out ${basename_id%%.*}Tss${disTSS}_body_Heatmap_360.pdf  --plotFileFormat pdf  --dpi 360  &
plotProfile -m matrix${disTSS}_${basename_id%%.*}_body.gz  -out ${basename_id%%.*}Tss${disTSS}_body_Profile.png   &
plotProfile -m matrix${disTSS}_${basename_id%%.*}_body.gz -out ${basename_id%%.*}Tss${disTSS}_Body_Profile_720.pdf --plotFileFormat pdf --perGroup --dpi 720  &
done
wait
