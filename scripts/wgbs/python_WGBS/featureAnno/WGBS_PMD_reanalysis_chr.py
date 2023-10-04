#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import math  # math.isnan()
import sys
import psutil
#from multiprocessing import Pool,Manager
import os
"""
WGBS results table
-> WGBS_PMD_reanalysis_chr  **now**
-> remodify 
-> R for geneId and genomic feature annotation

把基因组按照单位长度 segmentLength/stepInterval = 1000bp 切分每条染色体，并且结合甲基化测序文件，计算平均甲基化
必须分染色体进行，简化步骤不进行染色体核对，直接df.drop(索引)

segement each chromosome with 1000bp intervals, and calcu mean meathylation 

Note: WILL NOT CHECK THE CHROMOSOME ID!

"""
# def containNan(series):
#     return  any([ math.isnan(i) for i in series])

# #@func:
# #@input: start = ruleGenomeStart, stop = chromLen, step = stepInterval
# #@output: a list of tuple about each (start, end)
# def chromoInit(chromLen, stepInterval = 1000, ruleGenomeStart = 0):
#     nTimes,modBp = divmod(chromLen, stepInterval)
#     if modBp != 0:
#         start = [None] * (nTimes+1)
#         end = [None] * (nTimes+1)
#     else:
#         start = [None] * (nTimes)
#         end = [None] * (nTimes)
#     for index,i in enumerate(range(ruleGenomeStart, chromLen, stepInterval)):
#         start[index] = i
#         if index == nTimes:
#             end[index] = i+ modBp + 1 #右边取不到
#         else:
#             end[index] = i + stepInterval
#     return [ i for i in zip(start,end)]

# #@func:
# #@input: coordination list, sample name list for df construction
# #@output: a dataFrame with (start,end) and (samples)
# def construtDF_methySample(coordList,sampleNameList):
#     df = pd.DataFrame(columns=sampleNameList)
#     return pd.DataFrame(coordList, columns=["start","end"]).append(df)

# #@func: computation mean methylation
# #@input: (1).methylMatrix, (2).chromosome name, and (3).dataframe to be filled
# #@output: a dataFrame with meanMethylation 
# def meanMethyImpute_chr(chromosome, df_chr, methylMatrix):
#     for index,row in enumerate(zip(df_chr["start"], df_chr["end"])):
#         # match chr, start, end
#         equation = (methylMatrix["chr"] == chromosome ) & (methylMatrix["start"] >= row[0] ) & (methylMatrix["start"] < row[1])
#         if any(equation):
#             test = methylMatrix.loc[ equation ,  sampleList]
#             df_chr.loc[index,"A0316":"DH003"] = np.mean(test,axis=0)
#     return df_chr.loc[~df_chr.apply(containNan, axis=1) ]

# def meanMethyImpute_chr_muti(param):
#     index,row, chromosome, methylMatrix = param
#     dict_meanMeth = {}
#     # match chr, start, end
#     equation = (methylMatrix.chr == chromosome ) & (methylMatrix.start >= row[0] ) & (methylMatrix.start < row[1])
#     if any(equation):
#         test = methylMatrix.loc[ equation ,  sampleList]
#         dict_meanMeth[index] = np.mean(test,axis=0)
#         return dict_meanMeth

def meanMethyImpute_chr(methylMatrix,stepInterval):
    methylMatrix.loc[:,"segment"] = methylMatrix.start//stepInterval
    test=methylMatrix.groupby('segment').transform(np.mean)
    test.loc[:,"start"] = test.start//stepInterval*stepInterval
    test.loc[:,"stop"] = (test.stop//stepInterval+1)*stepInterval
    return test.drop_duplicates()

# #@func:  find non-density methylation region
# #@input: 
# #@output: 
# def findReduntInterval(secondMatrix, stepInterval = 1000, fold = 5,granularity_Step = 20):
#     down = []
#     upper = []
#     chromosome = []
#     num = 1     
#     for index in range(0, secondMatrix.shape[0], granularity_Step): #secondMatrix.shape[0]
#         startNew = secondMatrix.loc[int(index),"start"]
#         if ((startNew - num)//stepInterval >= fold) or (startNew/num > 5):
#             if index == secondMatrix.shape[0]  - 1:   # 最后一行
#                 break
#             else:
#                 neighbourDiff = np.diff(secondMatrix.loc[index: index+granularity_Step+1, "start"]) #差分 后-前
#             c = np.argmax(neighbourDiff)   # 按每行求出最大值的索引, c+1 - c 之间差距最大
#             # 计算 c/step+1 - c+1/step 之间的行索引
#             lowerIndex = secondMatrix.loc[ index+c, "start"] //stepInterval  +1  #表示后面一行开始
#             upperIndex = secondMatrix.loc[index+c+1, "start"]//stepInterval      #截止至前面一行
#             if upperIndex >= lowerIndex:
#                 down.append(lowerIndex)  # 索引是 down-1
#                 upper.append(upperIndex)   # 由于右侧取不到，所以直接使用
#         num = startNew
#     return [range(i[0],i[1]) for i in zip(down, upper)] 

# #@func:  drop non-density methylation region
# #@input: 
# #@output: 
# def simplyfyMatrix(reduntMatrix, intervaList):
#     indexRow = []
#     for i in intervaList:
#         indexRow.extend(list(i))
#     return reduntMatrix.drop(indexRow)

#@func:  combined extra file and define specific region
#@input: 
#@output: 
def fileExtraRegionAnno(meanMethyDF, dfAnnoCooorInfo,newColumnName ,stepInterval = 1000):
    meanMethyDF.loc[:, newColumnName] = False
    for row in zip(dfAnnoCooorInfo["start"], dfAnnoCooorInfo["end"]):
        front = row[0]//stepInterval
        tailD,tailM = divmod(row[1],stepInterval)
        if tailM == 0:
            tailD += 1
        else:
            tailD += 2
        meanMethyDF.loc[front:tailD, newColumnName] = True #"All" #2 warnings:See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy

sampleList = ['A0316', 'A0351', 'AH003', 'B0316', 'B0351', 'BH003', 'C0316', 'C0351', 'CH003', 'D0316', 'D0351', 'DH003']

parameter1 = sys.argv[1] # chromo
parameter2 = sys.argv[2] # chromo length
parameter3 = sys.argv[3] # num: one of PMD UMD LMD
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

methylMatrix = pd.read_table("../bsmkcov_AllSample_merged.txt")
# pmd_Stages = pd.read_table("Rwgbs_PMDIntersectAllStage.txt" ) 
# pmd_StageSetdiff = pd.read_table("Rwgbs_PMDSetdiffAllStage.txt" ) 
# basicLMR_stages = pd.read_table("basic_LMR_intersect_position4Homer.txt") 
# basicUMR_stages = pd.read_table("basic_UMR_intersect_position4Homer.txt") 


methylMatrix_chr = methylMatrix.loc[methylMatrix.chr == chromName].reset_index(drop=True)   #read WGBS methylation Matrix and filter sepecific chr
# pmd_Stages_chr = pmd_Stages.loc[pmd_Stages.seqnames == chromName].reset_index(drop=True)
# pmd_StageSetdiff_chr = pmd_StageSetdiff.loc[pmd_StageSetdiff.seqnames == chromName].reset_index(drop=True)
# basicLMR_stages_chr = basicLMR_stages.loc[basicLMR_stages.seqnames == chromName].reset_index(drop=True)
# basicUMR_stages_chr = basicUMR_stages.loc[basicUMR_stages.seqnames == chromName].reset_index(drop=True)

# 1. construct
# df_chr = construtDF_methySample(coordList= chromoInit(chromLen= chromLen, stepInterval = stepInterval), 
#             sampleNameList= sampleList)
# 2. simplify computation region
# deprecated：lack of accuracy and only 2-5% improve
#indexNAN = findReduntInterval(methylMatrix_chr,fold=1, granularity_Step=9, stepInterval = stepInterval)
#df_chr = simplyfyMatrix(df_chr, indexNAN)
#print("simplifying >>> reduce {} rows<<< ".format(len(indexNAN)))
# computute
#pmdMean1K_chr = meanMethyImpute_chr(chromosome = chromName,
#            df_chr=df_chr,methylMatrix=methylMatrix_chr )

meanMethyImpute_chr(methylMatrix_chr, stepInterval).to_csv("%dbp/wgbs_PMDStageSpecific_%s.csv" % (stepInterval,chromName),index=False)

print ('monitoring >>>当前进程的内存使用：%.4f GB<<<' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )
# # 3. basiic PMD， sepecific PMD
# fileExtraRegionAnno(pmdMean1K_chr, dfAnnoCooorInfo = pmd_Stages_chr, 
#             newColumnName = "BasicPMD", stepInterval = stepInterval) 
# fileExtraRegionAnno(pmdMean1K_chr, dfAnnoCooorInfo = pmd_StageSetdiff_chr, 
#             newColumnName = "SpecificPMD", stepInterval = stepInterval) 
# fileExtraRegionAnno(pmdMean1K_chr, dfAnnoCooorInfo = basicLMR_stages_chr, 
#             newColumnName = "BasicLMR", stepInterval = stepInterval)
# fileExtraRegionAnno(pmdMean1K_chr, dfAnnoCooorInfo = basicUMR_stages_chr, 
#             newColumnName = "BasicUMR", stepInterval = stepInterval)
# 4. export
#pmdMean1K_chr.loc[:,"start":"BasicUMR"].to_csv("%dbp/wgbs_PMDStageSpecific_%s.csv" % (stepInterval,chromName),index=False)
print("writing >>> wgbs_PMDStageSpecific_%s.csv <<< " % (chromName))

