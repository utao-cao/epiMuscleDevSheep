#!/bin/sh
###
 # @Authour: utao.cao
### 
#link: https://www.biostars.org/p/260702/
#%d %f
# 2015: SheepGenome_4.0.2015_MetaInfo.txt
# 2020: HuSheepGenomeMetaInfo.txt

# metaInfoGENOME can be downloaded from NCBI
# GenBank-Accn  RefSeq-Accn  Sequence-Length UCSC-style-name 

#if ($3-$2 >1) {
SOURCE_PATH="./wd_cyt/sheep/pj013_ATAC"
metaInfoGENOME=$1
echo $metaInfoGENOME
cd $SOURCE_PATH

# Function1: clip chromosome length
# modify
#cat  $metaInfoGENOME | awk '{printf("($1==\"%s\") {L=%d;B=int($2);E=int($3);S=$4;B=(B>=L?L-1:B);E=(E>=L?L-1:E);if($3 -$2 >1){printf(\"%s\\t%%d\\t%%d\\t%%f\\n\",B,E,S);next};}\n",$1,$3,$1);}' > ${SOURCE_PATH}/script.awk
# org 
#cat  $metaInfoGENOME | awk '{printf("($1==\"%s\") {L=%d;B=int($2);E=int($3);S=$4;B=(B>=L?L-2:B);E=(E>=L?L-1:E);if($3 - $2 >1){printf(\"%s\\t%%d\\t%%d\\t%%f\\n\",B,E,S);}next;}\n",$1,$3,$1);}' > ${SOURCE_PATH}/script.awk

# Function2: chromosome name among different versions
# 0. RefSeq <- GenBank
# 1: RefSeq -> GenBank
# -1. RefSeq -> chr-UCSC
# 2.  GenBank -> chr-UCSC
# example:
#RefSeq2GenBank=2
RefSeq2GenBank=3

# if 1st MEET the requirements( == RefSeq ); => ( GenBank-Accn) of metaInfoGENOME
if [ $RefSeq2GenBank == 1 ];then
    cat  $metaInfoGENOME |  awk '{printf("($1==\"%s\") {$1=\"%s\";printf($0\"\\n\");next;}\n",$2,$1);}' > ${SOURCE_PATH}/script_revertChromosome.awk
fi
# if 1st MEET the requirements( == GenBank-Accn );  =>  $4 ( chr-UCSC ) of metaInfoGENOME
if [ $RefSeq2GenBank == 2 ];then
    cat  $metaInfoGENOME |  awk '{printf("($1==\"%s\") {$1=\"%s\";printf($0\"\\n\");next;}\n",$1,$4);}' > ${SOURCE_PATH}/script_revertChromosome.awk
fi

# if 1st MEET the requirements( == GenBank-Accn );  =>  $4 ( chr-UCSC ) of metaInfoGENOME
if [ $RefSeq2GenBank == 3 ];then
    cat  $metaInfoGENOME |  awk '{printf("($1==\"%s\") {$1=\"%s\";printf($0\"\\n\");next;}\n",$4,$1);}' > ${SOURCE_PATH}/script_revertChromosome.awk
fi
