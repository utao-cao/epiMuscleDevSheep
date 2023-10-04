#!/bin/bash -e
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -J wgs
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
#https://www.jianshu.com/p/3a11c8f2498f
conda activate  ./conda_env/ngsGenome
CORE_FUL=64
CORE_LES=8
CORE_GRT=56
export GENOME_PATH="./reference/sheep/ncbi_sheep_2015"
SOURCE_DIR="./wd_cyt/sheep/sheep_WGS"
INPUT_DIR="${SOURCE_DIR}/fastq"
QC_DIR=${SOURCE_DIR}/output_qc
BWA_DIR=${SOURCE_DIR}/output_bwa
BAM_SOT=${BWA_DIR}/sorted
BAM_DUP=${BWA_DIR}/dedup

GVCF_DIR=${SOURCE_DIR}/output_gatk
VCF_DIR=${SOURCE_DIR}/output_gatk/population

#RGID=$3  ## Read Group，OR Lane ID
#library=$4  ## library ID
#sample=$5  ## sample ID
library='wgs'
RGID='00'

refGENOME="${GENOME_PATH}/GCA_000298735.2_Oar_v4.0_genomic.fa.gz"
BWA_INDEX=${SOURCE_DIR}/reference/bwadx_ovi4
SRR_meta=${INPUT_DIR}/SRR_metaInfo.txt
SRR_lis=${INPUT_DIR}/SRR_Acc_List.txt  

cd ${SOURCE_DIR}
DIRLISTs='fastq output_qc output_bwa/sorted output_bwa/dedup output_gatk/population reference '
for DIR_LIST in $DIRLISTs ;do
    if [ ! -d $DIR_LIST ];then
      mkdir -p $DIR_LIST
    fi
done

    if [  -f $refGENOME ];then
      pigz -d  $refGENOME
      
    fi
refGENOME="${GENOME_PATH}/GCA_000298735.2_Oar_v4.0_genomic.fa"
samtools dict  $refGENOME  -o $GENOME_PATH/GCA_000298735.2_Oar_v4.0_genomic.fa.dict

#step 01_1
#对下载的基因组建立索引，下面是建索引的脚本
#	-p STR   输出数据库的前缀；【默认和输入的文件名一致，输出的数据库在其输入文件所在的文件夹，并以该文件名为前缀。】
#	-a [is|bwtsw]   构建index的算法，有两个算法： is 是默认的算法，虽然相对较快，但是需要较大的内存，当构建的数据库大于
#               2GB的时候就不能正常工作了。 bwtsw 对于短的参考序列式不工作的，必须要大于等于10MB, 但能用于较大的基因组数据，比如人的全基因组。
bwa-mem2 index  -p $BWA_INDEX $refGENOME && echo "bwa index is done"  

#step 01_1
#https://www.biostars.org/p/147148/
#The SRA archive format ("vdb") contains an md5 checksum as well as a few other consistency checks (I think). 
#The sra-toolkit has a utility, vdb-validate which will report any errors in the data, 
#and perform an md5 checksum comparison.
#vdb-validate 
#fastq-dump --split-files 

#step 01_1
for sample in $(cat $SRR_lis);do
fastqc ${INPUT_DIR}/${sample}_1.fq.gz ${INPUT_DIR}/${sample}_2.fq.gz \
	-o ${QC_DIR}  && echo "${sample} fasqtc is done" > log_fastqc${i}.txt  &
done	

wait
#step 01_3
cd ${QC_DIR}
multiqc *.zip 


#step 02_1: align and sort
cat $SRR_meta |while read id;do 
	echo $id
	arr=($id)
	library=${arr[1]}
	sample=${arr[0]}
	#   -R $TAG \
  bwa-mem2 mem -t $CORE_FUL -M \
		$BWA_INDEX ${INPUT_DIR}/${sample}_1.fq.gz  ${INPUT_DIR}/${sample}_2.fq.gz  | \
		 ${SAMTOOLS} sort  -O bam  -@ ${CORE_FUL} -o - > $BAM_SOT/${sample}.sort.bam
done

#step 02_3
for sample in $( cat $SRR_lis  );do
##只是将重复标记，可以让后续软件识别，但是并没有删除
	gatk MarkDuplicates \
      -I $BAM_SOT/${sample}.sort.bam  \
      -M $BAM_DUP/${sample}.mk_metrics.txt \
      -O $BAM_DUP/${sample}.mk.bam && echo "mk ${sample} done" &
done
wait

#step 03_1
for sample in $( cat $SRR_lis );do
{
TAG="@RG\tID:$RGID\tPL:ILLUMINA\tSM:$sample\tLB:$library"

samtools addreplacerg -r $TAG  $BAM_DUP/${sample}.mk.bam -o $BAM_DUP/${sample}.mk_rg.bam
time samtools index  -@ ${CORE_FUL}  $BAM_DUP/${sample}.mk_rg.bam && echo "** ${sample}.mk.bam index done **" 

gatk HaplotypeCaller \
    --emit-ref-confidence GVCF \
    -R $refGENOME \
    -I $BAM_DUP/${sample}.mk_rg.bam  \
    --sample-name  ${sample}  \
    --native-pair-hmm-threads 8  \
     -O $GVCF_DIR/${sample}.HC.g.vcf.gz && echo "** GVCF ${sample}.HC.g.vcf.gz done **"

# 可以被删除清理的文件
#rm -f $BAM_SOT/${sample}.sort.bam 
} &

done
wait

#step 04: merge gvcf
# 按照","，把所有样本ID拆分出来存至数组中  $(echo $samples | tr "," "\n")
# 基于群体数据进行Joint genotyping, 同样有两种方式
## 第一，先合并所有的全gvcf结果，然后再统一进行GenotypeGVCFs，由于是全基因组，速度较慢
sample_gvcfs=""
outname='OLine6'
for sample in $(ls $GVCF_DIR/*.g.vcf.gz) ; do 
	sample_gvcfs=${sample_gvcfs}"-V ${sample}.HC.g.vcf.gz \\"\n
done

time gatk CombineGVCFs \
	-R $refGENOME \
	${sample_gvcfs} \
	-O $VCF_DIR/${outname}.HC.g.vcf.gz && echo "** ${outname}.HC.g.vcf.gz done ** " && \
time gatk GenotypeGVCFs \
	-R $refGENOME \
	-V $VCF_DIR/${outname}.HC.g.vcf.gz \
	-O $VCF_DIR/${outname}.HC.vcf.gz && echo "** ${outname}.HC.vcf.gz done ** "

#step 03_1 hard filter
#select SNP and perform filter
# zcat OLine6.HC.vcf.gz | grep -v "#"  | awk '{if($6 > 30) print $0}' | wc -l
# zcat OLine6.HC.vcf.gz | grep -v "#"  | awk -F\; '{match($0, "MQ=([.0-9]+);", k); if(k[1]> 30) print $0}'  | wc -l # qual/depth
gatk SelectVariants -select-type SNP -V $VCF_DIR/${outname}.HC.vcf.gz  -O $VCF_DIR/${outname}.HC.snp.vcf
gatk VariantFiltration -V $VCF_DIR/${outname}.HC.snp.vcf --filter-expression "QUAL < 30" --filter-name "PASS"  -O $VCF_DIR/${outname}.HC.snp.filter.vcf
