#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
file_abspath = '.'
def methylCor(dt,col_names):
    data = dt.loc[:,col_names].sort_index(axis=1)
    return data.corr()

dt = pd.read_csv('./test_BsmkCOV.txt',sep='\t',header=0)
corMat = methylCor(dt,col_names=dt.columns[3:])
print( corMat )
corMat.to_csv(os.path.join(file_abspath,"BsmkCOV_CpGCorrelation.txt"),index=False,sep='\t')
