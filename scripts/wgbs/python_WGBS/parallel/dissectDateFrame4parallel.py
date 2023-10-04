#!/usr/bin/env python 
# coding: utf-8

'''
objection: dissect data.frame

分割数据框
'''

import pandas as pd
import os
import argparse
# self-define
from createDir import create_dir
from pivotSummary_parallelVersion import reduceCategory,reduceDFmeno

#tmp_dir = '__TMP'
#df_Path = "pmdFullAnnoStat_100bp_4category.txt"

def get_parser():
    parser = argparse.ArgumentParser(description="若指定的文件夹不存在，则创建")
    parser.add_argument('df_path', metavar='DF_PATH', nargs = 1, type = str, help='指定文件')
    parser.add_argument('tmp_dir', metavar='TMP_DIR', nargs = 1, type = str, help='指定文件夹')
    parser.add_argument('n_cores', metavar='CORES', nargs = 1, type = int, help='核心数')
    return parser

def split(df_path, n_cores, tmp_dir):
    try:
        # read-in
        pmdAnnoStat = reduceDFmeno(pd.read_table(df_path,low_memory=False))
        cateCols = pmdAnnoStat.columns[12:]
        pmdAnnoStat.loc[:,cateCols] = reduceCategory(pmdAnnoStat, columns=cateCols)
        # pre-calcu
        nrows = pmdAnnoStat.shape[0]    #总行数
        averageRows = nrows//(n_cores)   
        lastCoresStartFrom = (nrows-averageRows*(n_cores-1))
        # segment
        segmentPoints = [i for i in range(0,nrows, averageRows)]
        intervaList = [ range(segmentPoints[index-1],segmentPoints[index]) for index in range(1,len(segmentPoints))] # 间隔 intervals
        segmentPoints[-1] = nrows
        print(nrows,averageRows,lastCoresStartFrom)
        print(intervaList)
        # output
        for index,item in enumerate(intervaList): 
            new_name = os.path.join(tmp_dir,'tmpDiseectedDF_{}'.format(index))
            pmdAnnoStat.iloc[item,:].to_csv(new_name,sep='\t',index=False)
    except Exception as e:
        print(e)

def main():
    args = vars(get_parser().parse_args())  # 返回对象的属性-属性值的字典对象  return a key-value dict
    df_path = args['df_path'][0]
    tmp_dir = args['tmp_dir'][0]
    n_cores = args['n_cores'][0]
    create_dir('.',tmp_dir)
    split(df_path, n_cores, tmp_dir)

if __name__ == '__main__':
    main()

