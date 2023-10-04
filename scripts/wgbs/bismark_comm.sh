#!/bin/sh -e
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -J bsmak1
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log

#####################################################
#                                                   #
#         INPUT INFORMATION REQUIRED                #
#                                                   #
#####################################################

#conda activate ~/conda_env/wgbs_wkflo/
export RUN_PATH=`pwd`
export SOURCE_DIR="./wd_cyt/ay_sheep/pj013_WGBS"
export WGBS_ENV_PATH="./conda_env/wgbs_wkflo/bin"
export GENOME_PATH="./reference/sheep/ncbi_sheep_2015" # 2015 

export CORE_NUM=4
export CORE_FUL=64
### Set the output dir to retain all the methylation called file:
fastqD="./wd_cyt/ay_sheep/pj013_WGBS/fastq"
bamD=${fastqD/fastq/align_bismark}
BAM_DIR=$bamD/bam_files  
UNSORTED_BAMdedupe_DIR=$bamD/unsort_Bismarkfile
METH_OUTPUT_DIR=$UNSORTED_BAMdedupe_DIR/methylation_extraction

#@rename  $pattern_str  $replace  $filename/dir  
rename_bysed(){
  pattern=$1
  replace=$2
  filename=$3
  if [ ! -f $filename ];then
    echo "file/dir does not exist"
  else {
    ls $filename  | grep --regexp=$pattern_str | sed 's@"$pattern"@"$replace"@g' | xargs -n 1 -I{} mv $filename {}
  }
  fi
}


if [ ! -d  $BAM_DIR ];then
  mkdir -p $BAM_DIR
fi

echo "prepare a dir contatining reference ending in .fa or .fasta(can be .gz) mandatory!"
#bismark_genome_preparation --genomic_composition   --parallel $CORE_FUL  $GENOME_PATH &
#wait

 	

JOB_NAME="bismark_genome_preparation"
 echo "FINISH:  $JOB_NAME"

#can set a var LIB storing list of IDs of files to be processed like ENCODE does; 
#and can mkdir for each sample a subfold under OUTPUT dir like bamD etc
#####################################################
#                                                   #
#         PIPELINE STARTS HERE                      #
#                                                   #
#####################################################
#rename files in propriate format and start from clean reads
	### Find the pair for each splittted fastq files to make it paired-end and feed it into the bismark alignment pipeline for further processing using trim_galore_bismark_alignment.sh script:
	### The fastq considered here would be of the format instance like : C6RLPANXX_s4_1_7bp_Index_7_SL83766.00 
	
	JOB_NAME="Rename files in propriate format: SampleID_trimTag_R1/R2.fq.gz"	
	for fastq in $( ls $fastqD/*fq.gz | grep  -v H003 ); do  #filenname with dir path
     
     sampleID=$( echo ${fastq##*/} | cut -f 1 -d _ | sed 's/��/(/' | sed 's/��/)/')     
     trimTag=$( echo ${fastq##*/} | cut -f 2 -d _)        #clean
     R1oR2=$( echo ${fastq##*/} | sed "s/.fq.gz//" | cut -f 3 -d _)   #1 or 2
     
             #echo   "${sampleID}_${trimTag}_${R1oR2}.fq.gz"
     #mv $fastq  ${fastq%/*}/${sampleID}_${trimTag}_${R1oR2}.fq.gz
     
 done
 
 JOB_NAME="Runnning bismark"
 
     sample_array01=(`ls ${fastqD} | grep "_clean_1.fq"  | grep  -v H003`)
     sample_array02=(`ls ${fastqD} | grep "_clean_2.fq"  | grep  -v H003`)
     
     arrayNumber=${#sample_array01[@]}
     step=$((arrayNumber/4+1))
  
  for i in $(seq -s' ' 0 $step);do
     
     echo "Start mapping turn:$i !"
     for j in $(seq 0 3);do
       
     	R1=${sample_array01[((i*4+j))]}
     	R2=${sample_array02[((i*4+j))]}
     	basename01=${R1##*/}
       basename02=${R2##*/} 
       echo   "basenames  is  ${basename01};"
       sampleID_R1=$( echo $basename01 | cut -f 1 -d _   )     
       sampleID_R2=$( echo $basename02 |cut -f 1 -d _  )     
       LIB=${sampleID_R1}
       trimTag=$( echo $basename01  | cut -f 2 -d _)        #clean
     
     	if [[ -f ${fastqD}/${R1} && -f ${fastqD}/${R2}  &&  ${sampleID_R1} = ${sampleID_R2} ]];then
           echo "${R1%_*} exist and $trimTag ;Will Start Mapping Soon: "
           echo "###################                     ${R1}  ${R2}"
           
           TEMP_DIR=${BAM_DIR}/${sampleID_R1}
           
           if [ -d $TEMP_DIR ];then
             rm -fr $TEMP_DIR
           else
             mkdir $TEMP_DIR
           fi
           
           echo "skip bismark"  #-p $CORE_NUM
           
           time $WGBS_ENV_PATH/bismark --bowtie2  --parallel 3 --non_directional --bam --temp_dir $TEMP_DIR --path_to_bowtie $WGBS_ENV_PATH/ \
                   -o ${BAM_DIR} $GENOME_PATH -1 ${fastqD}/$R1 -2 ${fastqD}/$R2   2>&1 >  ${BAM_DIR}/log_${sampleID_R1}.txt  &
			      #eval pid0${j}=$!
			     			      
           rm -rf ${TEMP_DIR}  && echo "sub  successfully"
			  fi
      done
      # ensure each parallel job finished !
      wait

   done
   wait 	
  	
  	
  JOB_NAME="Done deduplicate_bismark"
  
	cd $UNSORTED_BAMdedupe_DIR
	# if nums of samples < $CORE_FUL is ok
  for bam_file in $(/bin/ls $BAM_DIR/*_clean_1_bismark_bt2_pe.bam  | grep  -v H003 ); do
		#$WGBS_ENV_PATH/
  ( deduplicate_bismark -p --bam $bam_file  ) &
		sleep 20s
	done
  wait
  
  echo "FINISH:  $JOB_NAME"
  JOB_NAME="methylation_extractor"
  for bam_file in $(/bin/ls $UNSORTED_BAMdedupe_DIR/*_clean_1_bismark_bt2_pe.deduplicated.bam ); do
		#$WGBS_ENV_PATH/
		bismark_methylation_extractor -p --buffer_size 80% --merge_non_CpG --scaffolds --cytosine_report  -o ${METH_OUTPUT_DIR} --comprehensive --ucsc --gzip --parallel 20 \
                  --report --bedGraph --genome_folder $GENOME_PATH  $bam_file	&   # forget if proper to put in bg
    
  echo "FINISH:  $JOB_NAME"
	done
	
	wait