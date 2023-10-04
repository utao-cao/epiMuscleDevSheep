#!/bin/sh
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
# -*- coding: UTF-8 -*-
########################################################################
##
# @Description: Scripts to run transcript workflow
#
#@Author: CYT
#@Notion:
#1. need to rewrite  check_name_consistency()  mechanism to fit different names: ${read1%_*}
#@ref:
#1. seq使用变量，创建序列
#link: https://stackoverflow.com/questions/6191146/variables-in-bash-seq-replacement-1-10
#2. 动态创建变量： eval pid${j}=5
#link: https://www.cnblogs.com/bugutian/p/12030428.html
#3. while and 并列
#link: https://stackoverflow.com/questions/8239197/using-and-in-bash-while-loop
########################################################################
#Variable OR argument
sampleID=""
TotalSamples=27
cores=16  #4 part parallel

#directory : must end with /
fastqD=$1
#dir path + hisat2 index name;
indexGnomeD=$2
#fatsa file path
reference=$3
#gtf file path
annotationGtf=$4


bamD=${fastqD//fastq/bam}
cashD=${fastqD//fastq/cashAS}
countD=${fastqD//fastq/count}
#software
#hisat2 has been in base env
##samtools
export PATH="./software/samtools/bin/:$PATH"
########################################################################
##
# @Running
##
########################################################################
echo "create dir"
mkdir -p $fastqD $bamD $cashD $countD 
if [ ! -d ${indexGnomeD%/*} ];then
	mkdir ${indexGnomeD%/*}
	echo "need to create index"
	hisat2-build -p 64 $reference ${indexGnomeD} 
	fi
	
sample_array01=(`ls ${fastqD} | grep -v ^d | grep "_1.clean.fq"`)
sample_array02=(`ls ${fastqD} | grep -v ^d | grep "_2.clean.fq"`)

arrayNumber=${#sample_array01}
step=$((arrayNumber/4+1))
#reference: https://stackoverflow.com/questions/6191146/variables-in-bash-seq-replacement-1-10

# For RNA-seq reads, use "--dta/--downstream-transcriptome-assembly" 
# the command is similar to 'bowtie2'
#########################################
##
# @Mapping parallelly
###4 batch jobs to run parallel
########################################
for i in $(seq -s' ' 0 $step)
do
	echo "Start mapping turn:$i !"
	for j in $(seq 0 3)
	do
	read1=${sample_array01[((i*4+j))]}
	read2=${sample_array02[((i*4+j))]}
	if [ -f ${fastqD}${read1} ];then
        echo "${read1%_*} exist"
		if [ ${read1%_*} = ${read2%_*} ];then
			sampleID=${read1%_*}
			echo "Mapping ${read1%_*}  ${read1}  ${read2}"
		    time hisat2 -x ${indexGnomeD} -1 ${fastqD}${read1} -2 ${fastqD}${read2} -S ${bamD}${sampleID}.sam -p $cores --dta  2>&1 | tee ${bamD}${sampleID}.stat &
			eval pid0${j}=$!
		fi
	fi
	done
# ensure each parallel job finished !
	wait $pid00; wait $pid01; wait $pid02; wait $pid03
	unset $pid00; unset $pid01; unset $pid02; unset $pid03
	
    echo "Start with: "$(date "+%Y-%m-%d %H:%M:%S")
    #wait
    echo "sleep 30s"; sleep 30
    echo "end with: "$(date "+%Y-%m-%d %H:%M:%S")
	echo "wait Mapping down: $i"
	
#########################################
# @Sort And Count Reads
# also can be modidied to paralle
# All dir must names by suffix "D"
########################################
	for samFile in `ls $bamD | grep sam`
	do
		sampleID=${samFile%.*}
		time samtools sort -@ 64 -o ${bamD}${sampleID}.sorted.bam  ${bamD}${samFile}
		echo "delet ${samFile}"
		rm ${bamD}${samFile}
	
		wait ; sleep 10
	
		echo "Count ${sampleID}.sorted.bam"
		time featureCounts -T 64 -p -t exon -g gene_id -a $annotationGtf -o ${countD}${sampleID}_featureCounts.txt \
			${bamD}${sampleID}.sorted.bam 2>&1 ${countD}${sampleID}.summary.txt

		wait ; sleep 10
	done

	wait
done

echo "Congratulation:All Sample Mapping Finished, Sorted And Reads Counted!"

