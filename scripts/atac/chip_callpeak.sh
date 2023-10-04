#!/bin/sh
###
 # @Authour: utao.cao
### 
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
 mkdir sra fastq bam log 
 sample="SRR5070526"
 BT2_INDEX="../reference/Sheep_bowtie"
 fq1=

 GSIZE=$(bowtie2-inspect $BT2_INDEX |perl -ne 'BEGIN{$n=0} next if(/^>/);s/[^ATGC]//gi; $n+=length($_); END{print int($n*0.85);}')  &


 fastp -o ./fastq/clean/${sample}_1.fastq.gz -i ./fastq/${sample}_1.fastq \
#        -O ./fastq/clean/${sample}_2.fastq.gz -I ./fastq/${sample}_2.fastq \
        -Q --thread=16 --length_required=50 --n_base_limit=6 --compression=6 &

cd ~/seqData


 pid=$!
 wait $pid
 echo " Completed Index !\n"
 time bowtie2  -p ${CORE_FUL}  --very-sensitive -X 2000 -x  $BT2_INDEX -U ./fastq/clean/${sample}_1.fastq.gz  | \
          samtools view -bhS -q 30 | \
          ${SAMTOOLS} sort -l 9 -O bam  -@ ${CORE_FUL} -o - > ./bam/${sample}.highQuality_sort.bam
 echo "mapping and sort\n"
 
 time  ${SAMTOOLS} index -@ ${CORE_FUL} ./bam/${sample}.highQuality_sort.bam
 
## option 1 
### firstly we just keep the high mapping quality reads according to ENCODE project guideline.
#samtools view -bhS -q 30  ${id%%.*}.sam > ${id%%.*}.highQuality.bam  
### -F 1548 https://broadinstitute.github.io/picard/explain-flags.html 
#samtools sort   ${id%%.*}.highQuality.bam ${id%%.*}.highQuality.sorted  ## prefix for the output   
#samtools index ${id%%.*}.highQuality.sorted.bam 
## option 2
# ## Then we just keep the unique mapping reads according to the majority tutorial.
# grep -v "XS:i:" ${id%%.*}.sam |samtools view -bhS - >${id%%.*}.unique.bam
# samtools sort   ${id%%.*}.unique.bam ${id%%.*}.unique.sorted  ## prefix for the output   
# samtools index ${id%%.*}.unique.sorted.bam 
 wait
 
 echo "Effective Genome Size" $GSIZE
 macs2 callpeak -t ./bam/${sample}.highQuality_sort.bam  --keep-dup all --nomodel -- extsize 100 -f BAM \
        -g $GSIZE -n ${sample} 2> ./log/${sample}.masc2.log
