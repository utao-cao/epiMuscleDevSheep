#!/bin/sh -e
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -J valSplit
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
#cat SRR_Acc_List.txt | while read id; do 
#  echo $id;
#  nohup sh -c  " prefetch $id --max-size 200G  && \
#      vdb-validate  $sample  && \
#      fastq-dump --split-files ${sample}/${sample}.sra  " > ${id}.txt & ;
#done
#wait

for sample in $(cat SRR_Acc_List.txt);do
  {
  echo $sample
    vdb-validate  $sample && \
  		fasterq-dump --threads 16 --split-3  -p   ${sample}/${sample}.sra
  } &
done
#
#cat SRR_metaInfo.txt  |while read id;do 
#  echo $id
#  arr=($id)
#  library=${arr[1]}
#  sample=${arr[0]}
#  
#  echo $library  echo $sample
#  fasterq-dump --threads 58 --split-3  -p   ${sample}/${sample}.sra  
#  done

pigz -p 58  *fastq

for i in $(ls *fastq.gz);do mv $i  $(echo ${i} | sed 's/.sra//g' | sed 's/.fastq.gz/.fq.gz/g') ;done

wait
echo "done"
#SRR10821772
#SRR11657579
#SRR11657580
#SRR11657581
#SRR11657582
#SRR11657583
#SRR11657683
#SRR11657684
#SRR11657685
#SRR11657686
#SRR11657687