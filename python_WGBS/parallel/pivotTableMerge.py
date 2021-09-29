#!/usr/bin/env python 
# coding: utf-8
import numpy as np
import pandas as pd
import sys
import os
import re
import argparse
from functools import reduce
from createDir import create_dir
from pivotSummary_parallelVersion import reduceCategory,reduceDFmeno,mem_usg

pd.set_option('precision', 4)
def get_parser():
    parser = argparse.ArgumentParser(description="若指定的文件夹不存在，则创建")
    parser.add_argument('df_path', metavar='DF_PATH', nargs = 1, type = str, help='指定文件路径')
    parser.add_argument('out_name', metavar='OUT_NAME', nargs = '+', type = str, help='指定输出文件')
    return parser

def find_byPattern(df_path,pattern): 
    if os.path.exists(df_path):
        #print(os.listdir(df_path))
        namePatn = re.compile(pattern)
        f_nameList = [] 
        for root,dirs,file in os.walk(df_path):
           f_nameList.extend( [ os.path.join(root,f) for f in file if namePatn.match(f)] )
        return f_nameList
        
def readTableList(df_paths):
    df_list = []
    for i in df_paths:
        df = reduceDFmeno( pd.read_table(i, low_memory=False) )
        cateCols = df.columns[12:]
        df.loc[:,cateCols] = reduceCategory(df, columns=cateCols)
        df_list.append(df)
    return df_list

def main():
    pattern = r"tmp[a-zA-Z]*_[0-9]*meanMethydescribe.txt$"  # _individul
    args = vars(get_parser().parse_args())  # 返回对象的属性-属性值的字典对象  return a key-value dict
    df_path = args['df_path'][0]
    out_name = args['out_name'][0]
    df_list = readTableList(find_byPattern(df_path,pattern))
    mergeDF = pd.concat(df_list)
    mergeDF.loc[mergeDF.duplicated(subset=['geneId','group','cateRange','typeMethy'], keep=False)].to_csv("{}Duplicatedmerge.txt".format(out_name),sep='\t')
    mergeDF.to_csv("{}merge.txt".format(out_name),sep='\t')

if __name__ == '__main__':
    main()