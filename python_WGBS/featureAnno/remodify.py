#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import math  # math.isnan()
import sys
import os
import re
"""
WGBS results table
-> WGBS_PMD_reanalysis_chr
-> remodify **now**
-> R for geneId and genomic feature annotation

1. 把基因组按照单位长度 segmentLength/stepInterval = 1000bp 直接根据坐标运算，取余（向量化处理），得到segment
2. 根据提供的注释区间文件，进行相对应注释：PMD，UMD,LMR,差异DMR等
3. 改进速度

技术要点：
1. zip(df.to_numpy().tolist()).__next__()
2. pandas 遍历提速iloc 太慢了：iat/at, iterrows, apply, 列表构造， 数pandas or numpy组操作 df['a'] or df['a'].values()


"""
# parameter
parameter1 = sys.argv[1] # chromo
parameter2 = sys.argv[2] # chromo length
parameter3 = sys.argv[3] # stepInterval

# regex
groupRegex = re.compile(r'[A-D]{2}')
diffTypeRegex = re.compile(r'DM[RL]?')

def findOverlap(start1,start2,end1, end2):
    return max(start1,start2) < min(end1, end2)


'''
只适合两者都是顺序的
第一步：approximate 对于每一个(start-End)，找到起点，计算大致区间
第二步：iterate 迭代区间:  区间连续，从起点开始-
第三步：start new turn from the terminator of last iteration  
      返回 2.精确结束区间作为起点（LastIndex是真实索引, 不被改变）

否则，需要全部迭代；
'''
def boostFindOverlap_2steps(targetStartEnd, meanMethyDF, nStart, stepInterval):
    Slist = []
    count=0
    # as for stage-specific region, nStart >> nTerm  [meanMethyDF not fully ordered]
    for FIndex,FrowMethy in enumerate(zip(meanMethyDF.loc[nStart:,"start"], meanMethyDF.loc[nStart:,"stop"])): # 1st time, evaluate the approximal range
        if findOverlap(FrowMethy[0],targetStartEnd[0],FrowMethy[1],targetStartEnd[1]):
            count += 1
            nStart = FIndex
            nTerm = targetStartEnd[1]//stepInterval + 2  # treat as the maximum condition
            if nStart >= nTerm:      #more conservative way
                nTerm = nStart + (targetStartEnd[1] - targetStartEnd[0])//stepInterval + 2             
            #print(targetStartEnd,nStart,nTerm) 
            break #failed      # once get approximal range, break in time
    if count == 0:
        print("not found{}".format(targetStartEnd)) # 一直没找到, return None
    else:
        for SIndex,SrowMethy in enumerate(zip(meanMethyDF.loc[nStart:nTerm,"start"], meanMethyDF.loc[nStart:nTerm,"stop"])): # 2nd time, iteration among the interval
            if findOverlap(SrowMethy[0],targetStartEnd[0],SrowMethy[1],targetStartEnd[1]):
                Slist.append(FIndex + SIndex + 1)     # index calculate
        return range(nStart, Slist[-1])

def locIntervalPoint(startB, endB, stepInterval): 
    front = startB//stepInterval*stepInterval
    tailD,tailM = divmod(endB,stepInterval)
    if tailM == 0:
        tailD += 2   #是否应该算？ 还需要考虑仔细
    else:
        tailD += 1
    return (front,tailD*stepInterval)  # return coordination

#@func:  combined extra file and define specific region
#@input: 
#@output: 
def fileExtraRegionAnno(meanMethyDF, annoCoorInfoDF,newColname = '',stepInterval = 1000, annoInnerColname=False, bothOrder = False):
    #meanMethyDF.loc[:, newColname] = False
    indexRangeList = [[0]]
    if isinstance(annoInnerColname,(str,list)):  # argue2: tuple wrapped - OR
        for row in zip(annoCoorInfoDF["start"], annoCoorInfoDF["end"], annoCoorInfoDF[annoInnerColname].to_numpy().tolist()):   # load extra columns info
            startEnd = locIntervalPoint(row[0],row[1],stepInterval)   # calculate start and end segments
            #print(startEnd)
            if bothOrder:      #结束区间作为起点
                indexRangeList.append(boostFindOverlap_2steps(startEnd, meanMethyDF, indexRangeList[-1][-1], stepInterval))
            else:       # 常规只取一段大概区间进行迭代
                indexRangeList.append(boostFindOverlap_2steps(startEnd, meanMethyDF, 0, stepInterval))
            if not indexRangeList[-1]:
                indexRangeList.pop()
            else:
                meanMethyDF.loc[indexRangeList[-1], newColname] = row[2] # allow rename column; start from 3rd iterm of tuple 
    else:  # PMD,UMD,LMD like
        for row in zip(annoCoorInfoDF["start"], annoCoorInfoDF["end"]): # 列表构造
            startEnd = locIntervalPoint(row[0],row[1],stepInterval)
            if bothOrder:      #结束区间作为起点
                print(indexRangeList)
                indexRangeList.append(boostFindOverlap_2steps(startEnd, meanMethyDF, indexRangeList[-1][-1], stepInterval))
            else:       # 常规只取一段大概区间进行迭代
                indexRangeList.append(boostFindOverlap_2steps(startEnd, meanMethyDF, 0, stepInterval))
            if not indexRangeList[-1]: # 能够在meanMethy上找到的情况下
                indexRangeList.pop()
            else:
                meanMethyDF.loc[indexRangeList[-1], newColname] = True #"All" #2 warnings:See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy

def fixChrColumn(df):
    df.chr = ['chr{}'.format(i) for i in df.chr]

sampleList = ['A0316', 'A0351', 'AH003', 'B0316', 'B0351', 'BH003', 'C0316', 'C0351', 'CH003', 'D0316', 'D0351', 'DH003']
dmrList = ["dmrs_callDMRBSobj_AB.txt", "dmrs_callDMRBSobj_AC.txt", "dmrs_callDMRBSobj_AD.txt", 
           "dmrs_callDMRBSobj_BC.txt", "dmrs_callDMRBSobj_BD.txt", "dmrs_callDMRBSobj_CD.txt"]

if isinstance(parameter1, str):
    chromName = parameter1.strip()
else:
    raise Exception('{} 作为第一个参数：染色体号<str>；仅支持UCSC style<chr1> '.format(parameter1))

if isinstance(parameter2, str):
    chromLen = int(parameter2)
else:
    raise Exception('{} 作为第二个参数：染色体长度<int/str>'.format(parameter2)) 

stepInterval = int(parameter3)
print('calling >>> python WGBS_PMD_reanalysis_chr.py %s %d %d <<<' % (chromName, chromLen, stepInterval))


pmd_Stages = pd.read_table("Rwgbs_PMDIntersectAllStage.txt" ) 
pmd_StageSetdiff = pd.read_table("Rwgbs_PMDSetdiffAllStage.txt" ) 
basicLMR_stages = pd.read_table("basic_LMR_intersect_position4Homer.txt") 
basicUMR_stages = pd.read_table("basic_UMR_intersect_position4Homer.txt")

pmdA = pd.read_table("data/Rwgbs_pmdA.txt")
pmdB = pd.read_table("data/Rwgbs_pmdB.txt")
pmdC = pd.read_table("data/Rwgbs_pmdC.txt")
pmdD = pd.read_table("data/Rwgbs_pmdD.txt")

lmrA = pd.read_table("data/Rwgbs_lmrA.txt")
lmrB = pd.read_table("data/Rwgbs_lmrB.txt")
lmrC = pd.read_table("data/Rwgbs_lmrC.txt")
lmrD = pd.read_table("data/Rwgbs_lmrD.txt")

umrA = pd.read_table("data/Rwgbs_umrA.txt")
umrB = pd.read_table("data/Rwgbs_umrB.txt")
umrC = pd.read_table("data/Rwgbs_umrC.txt")
umrD = pd.read_table("data/Rwgbs_umrD.txt")

pmdMean1K_chr = pd.read_csv("{}bp/wgbs_PMDStageSpecific_{}.csv".format(stepInterval,chromName))   #read WGBS methylation Matrix and filter sepecific chr
pmd_Stages_chr = pmd_Stages.loc[pmd_Stages.seqnames == chromName].reset_index(drop=True)
pmd_StageSetdiff_chr = pmd_StageSetdiff.loc[pmd_StageSetdiff.seqnames == chromName].reset_index(drop=True)
basicLMR_stages_chr = basicLMR_stages.loc[basicLMR_stages.seqnames == chromName].reset_index(drop=True)
basicUMR_stages_chr = basicUMR_stages.loc[basicUMR_stages.seqnames == chromName].reset_index(drop=True)

pmdA_chr = pmdA.loc[pmdA.seqnames == chromName].reset_index(drop=True)
pmdB_chr = pmdB.loc[pmdB.seqnames == chromName].reset_index(drop=True)
pmdC_chr = pmdC.loc[pmdC.seqnames == chromName].reset_index(drop=True)
pmdD_chr = pmdD.loc[pmdD.seqnames == chromName].reset_index(drop=True)

lmrA_chr = lmrA.loc[lmrA.seqnames == chromName].reset_index(drop=True)
lmrB_chr = lmrB.loc[lmrB.seqnames == chromName].reset_index(drop=True)
lmrC_chr = lmrC.loc[lmrC.seqnames == chromName].reset_index(drop=True)
lmrD_chr = lmrD.loc[lmrD.seqnames == chromName].reset_index(drop=True)

umrA_chr = umrA.loc[umrA.seqnames == chromName].reset_index(drop=True)
umrB_chr = umrB.loc[umrB.seqnames == chromName].reset_index(drop=True)
umrC_chr = umrC.loc[umrC.seqnames == chromName].reset_index(drop=True)
umrD_chr = umrD.loc[umrD.seqnames == chromName].reset_index(drop=True)

#
# 3. basiic PMD， sepecific PMD
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = pmdA_chr, 
            newColname = "A_PMD", stepInterval = stepInterval) 
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = pmdB_chr, 
            newColname = "B_PMD", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = pmdC_chr, 
            newColname = "C_PMD", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = pmdD_chr, 
            newColname = "D_PMD", stepInterval = stepInterval) 
            
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = lmrA_chr, 
            newColname = "A_LMR", stepInterval = stepInterval) 
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = lmrB_chr, 
            newColname = "B_LMR", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = lmrC_chr, 
            newColname = "C_LMR", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = lmrD_chr, 
            newColname = "D_LMR", stepInterval = stepInterval) 
            
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = umrA_chr, 
            newColname = "A_UMR", stepInterval = stepInterval) 
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = umrB_chr, 
            newColname = "B_UMR", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = umrC_chr, 
            newColname = "C_UMR", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = umrD_chr, 
            newColname = "D_UMR", stepInterval = stepInterval) 
            
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = pmd_Stages_chr, 
            newColname = "BasicPMD", stepInterval = stepInterval) 
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = basicLMR_stages_chr, 
            newColname = "BasicLMR", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = basicUMR_stages_chr, 
            newColname = "BasicUMR", stepInterval = stepInterval)
fileExtraRegionAnno(pmdMean1K_chr, annoCoorInfoDF = pmd_StageSetdiff_chr, 
            newColname = "SpecificPMD", stepInterval = stepInterval) 

for i in dmrList:
    sample = groupRegex.findall(i)[0]
    diffType = diffTypeRegex.findall(i)
    df = pd.read_csv(i)
    fixChrColumn(df)
    fileExtraRegionAnno(pmdMean1K_chr, 
            annoCoorInfoDF = df.loc[df.chr == chromName].reset_index(drop=True), annoInnerColname=["areaStat","diff.Methy"],
            newColname = ["{}_{}".format(sample,x) for x in ["areaStat","diff.Methy"]] , stepInterval = stepInterval)

# 4. export
path="remodi{}bp".format(stepInterval)
if not os.path.exists(path):
    os.mkdir(path)
pmdMean1K_chr.to_csv("remodi%dbp/wgbs_PMDStageSpecific_%s.csv" % (stepInterval,chromName),index=False)
print("writing >>> wgbs_PMDStageSpecific_%s.csv <<< " % (chromName))