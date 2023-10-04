#!bin/python
#-*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
from functools import reduce

'''
merge multisamples and calculate cor base on methylation of CpG sites

additional operation:
1. drop zeros in one row after merging all the cov files through --gapThreshold

Memo: can be optimised to reduce memory occupation like ./parallel/pivotalSummary...

'''
def parse_args():
    parser = argparse.ArgumentParser(prog="merge multisamples and calculate cor base on CpG sites",prefix_chars="-",
        description="the message into before help: 2021/03/29",
        epilog='coverage files sep by \t, with colnames: chr     start   stop    sampleName')
    parser.add_argument('-o','--outdir',dest='outdir')
    parser.add_argument('-f','--files',dest='coverage_file', type=str,nargs='+',help='list')
    parser.add_argument('-g','--gapThreshold',dest='gapThreshold', type=float)
    return parser.parse_args()
    
# @fun: cal pearson's correlation 
# @return: correlation matrix
# @parameter:
# default method='pearson' | alternatively "spearman"
def calculate_methyl_cor(dt,col_names):
    data = dt.loc[:,col_names].sort_index(axis=1)
    return data.corr()

# @fun: counter num of `targetItem` (designed for row iteration) and which row to drop`` 
# @return: unqualified row index
# @parameter:   
# `gapRateLessThan` range(0,1): retain threshold
# # 1 == 100%  only when gapped rate%(`targetItem`) >= 100% , drop  
# #    =>>> keep all rows     
# # 0.3 >= 0   only when gapped rate% >=  0.3, drop
# #    =>>> keep rows which gapped rate% <  0.3
def zeroDropper_rowise(index, row, targetItem, dataColumns, gapRateLessThan=1, exclude_columns = 3):
    counter = 0
    for j in row[exclude_columns-1:]:
        if j == targetItem: 
            counter +=1
    if counter/dataColumns >= gapRateLessThan:
        return index

# @fun:  a counter: recoder the maximum consecutive `targetItem`, 
# @return:  a <dict> {row:maximun count}
# @parameter:   
def consecZeroDropper_rowise(index, row, targetItem, dataColumns, exclude_columns = 3):
    consec_counter = 0
    lastItem = 0
    consecutive_dict={}
    for j in row[exclude_columns-1:]:
        if j == targetItem and lastItem != targetItem:
            consec_counter = 1
        elif j == targetItem and lastItem == targetItem:
            consec_counter += 1      # update a new counter
        elif j != targetItem and lastItem == targetItem:
            consecutive_dict[index] = consec_counter
            if consec_counter and max(consecutive_dict.values())/dataColumns >= 0.3:
                pass #print( "%s\t%d\t%d" % (index, row[1], max(consecutive_dict.values())) )
        lastItem = j
    return consecutive_dict
    
# @fun:  wrapped func to rowwisely drop lines whose gap_rate more than {gap_rate_threshold}%
# @return:  2 objs: unqualified row index; consecutive counter results<dict>
# @parameter:   
# targetItem: 0   0 coverage 
# exclude_columns: "seqnames start end" these columns are not for counting
def zeroDropper(df_num, targetItem, exclude_columns = 3, gapRateLessThan = 1):
    dropLines=[]
    consecutive_dict={}
    columns = df_num.shape[1] - exclude_columns
    for index,row in enumerate(df_num.values):
        dropIndex = zeroDropper_rowise(index, row, targetItem, dataColumns=columns, gapRateLessThan = 1, exclude_columns = 3)
        consectDict = consecZeroDropper_rowise( index, row, targetItem, dataColumns=columns, exclude_columns = 3)
        if dropIndex:
            dropLines.append( dropIndex )
        if any(consectDict.values()):    
            consecutive_dict.update( consectDict )
    
    return (dropLines,consecutive_dict)

# @fun:  reduce + merge
# @return:  merged 
def reduce_merge(file_abspath, file_list): 
    df_list = []
    for i in file_list:
        df = pd.read_csv(os.path.join(file_abspath,i),sep="\t") #skip first line 跳过第一行, skiprows=1,nrows=5000000
        df.iloc[:,3:]=df.iloc[:,3:].round(2)
        df_list.append(df)
    df_merge = reduce(lambda left, right: pd.merge(left, right,on = ['chr', 'start', 'stop'],how =  'outer', sort=False), df_list).fillna(0)
    return df_merge

def main():
    args = parse_args()
    output_file = args.outdir
    countL = args.coverage_file    # a list of files
    gapRateThreshold = args.gapThreshold
    file_relaPath='.'    
    file_abspath=os.path.abspath(file_relaPath) #get abs dir path 获取文件的绝对路径
    #dir_path = os.path.dirname(file_abspath)
    #list out files, and assign suffix 列出文件，并指定后缀名
    #countL = [x for x in os.listdir(file_abspath) 
    # if os.path.isfile(os.path.join(file_abspath,x)) and os.path.splitext(x)[1]=='.bed']
    print(countL)

    df_merge = reduce_merge(file_abspath, countL)
    #print(df_merge.head(5))
    dropResults = zeroDropper(df_merge, targetItem = 0, gapRateLessThan= gapRateThreshold)
    consectutiveDict = dropResults[1]
    dropLines = dropResults[0]
    if dropLines:
        df_final = df_merge.drop( dropLines , axis=0)
    else:
        #print(dropLines)
        df_final = df_merge 
    # dir_path:  output to .. 
    # file_abspath: output to .
    # not output column names and row index
    df_final.to_csv(os.path.join(file_abspath, "%s_merged.txt"  % output_file), index=False,sep='\t')
    df_consectutive = pd.DataFrame.from_dict(consectutiveDict, orient ='index').reset_index().rename( columns = {'index':'row'})
    df_consectutive.to_csv(os.path.join(file_abspath,"%s_tandemedRows.txt"  % output_file ), index=False, sep='\t')


    # sample correlation
    corMat = calculate_methyl_cor(df_final,col_names=df_final.columns[3:])
    corMat.to_csv(os.path.join(file_abspath,"%s_CpGCorrelation.txt" % output_file ),index=False,sep='\t')
    print( corMat )

    #print(df_final.columns[3:])

if __name__ == "__main__":
    main()