#!/usr/bin/env python 
# coding: utf-8
import numpy as np
import pandas as pd
import math  # math.isnan()
import sys
import os
import re
import argparse
from createDir import create_dir

pd.set_option('precision', 4)

'''
we aim to calculate mean methylation of different methlyType for each gene, which looks like:
gene  group type meanMethylation
gene1  A    pmd   60.8
gene1  A    umr   30.5
gene1  B    lmr   40.0

regular workflow:
1. read in wide-form dataframe(pre-segmented files) optimising the memory usage(@reduceDFmeno);
2. flatten the methylation type(XPMD XLMR  XUMR) and  @keepConsistentGroup;
3. #comm exclude ambiguous annotation -> insufficient genes to build model
4. define high methylated region(HMD) and ambiguous merhylation region(AMD) (updated: 2021-09-28);
5. flatten mean Methy records for every group(A B C D) and  @keepConsistentGroup;
6. melt all the methylation types(5) and groupby-summary-export

Exceptional circumstances:
but there are AMDs which annotated by more than 2 types( we defined these region alone);

coding issues:
1. don't forget to reset_index before iteration
2. np.nan != np.nan
'''

def get_parser():
    parser = argparse.ArgumentParser(description="若指定的文件夹不存在，则创建")
    parser.add_argument('df_path', metavar='DF_PATH', nargs = 1, type = str, help='指定文件')
    parser.add_argument('out_dir', metavar='OUT_DIR', nargs = '+', type = str, help='指定输出文件夹')
    return parser

# keep two columns in one df consistent(only keep col1) 
def keepConsistentGroup(DF,col1,col2):
    groupDropList= []
    for index,row in enumerate(zip(DF[col1],DF[col2])):
        if len(set(row)) > 1:
            groupDropList.append(index)
    return DF.drop(groupDropList, axis=0).drop(col2, axis=1)

def reduceCategory(largeDF, columns):
    # 检查是否包含
    if len(columns) == 0:
        columns = largeDF.columns
    convert_obj = pd.DataFrame()
    for col in columns:
        num_uniqVal = len(largeDF[col].unique())
        num_totalVal = len(largeDF[col])
        if num_uniqVal/num_totalVal < 0.5:
            convert_obj.loc[:,col] = largeDF[col].astype('category')
        else:
            convert_obj.loc[:,col] = largeDF[col]
    return convert_obj

def reduceDFmeno(largeDF):
    print(mem_usg(largeDF))
    # int
    converted_int = largeDF.select_dtypes(include=['int']).apply(pd.to_numeric, downcast='unsigned')
    largeDF[converted_int.columns] = converted_int
    # float
    converted_float = largeDF.select_dtypes(include=['float']).apply(pd.to_numeric, downcast='float')
    largeDF[converted_float.columns] = converted_float

    print(mem_usg(largeDF))
    return largeDF

def containNan(series, how='any'):
    res = []
    for i in series:
        if how == 'any':
            if math.isnan(i):
                res = True
                break
            else:
                res.append(math.isnan(i))
        if how == 'all':
            if not math.isnan(i):
                res = False
                break
            else:
                res.append(math.isnan(i))
    if isinstance(res, bool) == 1:
        return res
    elif how == 'any':
        return any(res)
    elif how == 'all':
        return all(res)

def mem_usg(pdobj):
    if isinstance(pdobj, pd.DataFrame):
        usage_b = pdobj.memory_usage(deep=True).sum()
    else:
        usage_b = pdobj.memory_usage(deep=True)
    usage_mb = usage_b/1024**2
    return "{:03.2f} MB".format(usage_mb)

def main():
    args = vars(get_parser().parse_args())  # 返回对象的属性-属性值的字典对象  return a key-value dict
    df_path = args['df_path'][0]  #'parallel/tmpDiseectedDF_1025' #
    out_dir = args['out_dir'][0]  #'parallel/output' #
    create_dir('.',out_dir)

    pmdAnnoStat = reduceDFmeno( pd.read_table(df_path, low_memory=False) )
    cateCols = pmdAnnoStat.columns[12:]
    pmdAnnoStat.loc[:,cateCols] = reduceCategory(pmdAnnoStat, columns=cateCols)
    print(mem_usg(pmdAnnoStat))

    # category columns names by diff annotations
    sampleName = ('A0316','A0351','AH003','B0316','B0351','BH003','C0316','C0351','CH003','D0316','D0351','DH003')
    diffMethy = set([ re.match("[A-D]*_diff.Methy",x)[0]  for x in pmdAnnoStat.columns if re.match("[A-D]*_diff.Methy",x) ])
    pmdType = set([ re.match("^[BS][a-z]+[PMDURL]+",x)[0]  for x in pmdAnnoStat.columns if re.match("^[BS][a-z]+[PMDURL]+",x) ])
    areaStat = set([ re.match("[A-D]*_areaStat",x)[0]  for x in pmdAnnoStat.columns if re.match("[A-D]*_areaStat",x) ])
    pmd = set([ re.match("[A-D]_PMD",x)[0]  for x in pmdAnnoStat.columns if re.match("[A-D]_PMD",x) ])
    lmr = set([ re.match("[A-D]_LMR",x)[0]  for x in pmdAnnoStat.columns if re.match("[A-D]_LMR",x) ])
    umr = set([ re.match("[A-D]_UMR",x)[0]  for x in pmdAnnoStat.columns if re.match("[A-D]_UMR",x) ])

    dropCols=list(areaStat.union(diffMethy).union(pmdType))
    pmdAnnoStatSimp = pmdAnnoStat.drop(dropCols,axis=1, inplace=False)

    # flatten the PMD LMR  UMR
    ## wide form to longer form: all columns but pmd-related keep still, melts pmd
    pmdAnnoStatSimp = pmdAnnoStatSimp.melt(id_vars=set(pmdAnnoStatSimp.columns).difference(pmd), 
                                var_name="pmdGroup", value_name = 'pmd')
    pmdAnnoStatSimp.loc[:,'pmdGroup'] = pmdAnnoStatSimp['pmdGroup'].str[0:1]
    ## the same treatments apply to lmr
    pmdAnnoStatSimp = pmdAnnoStatSimp.melt(id_vars=set(pmdAnnoStatSimp.columns).difference(lmr), 
                                var_name="lmrGroup", value_name = 'lmr')
    pmdAnnoStatSimp.loc[:,'lmrGroup'] = pmdAnnoStatSimp['lmrGroup'].str[0:1]
    ## drop inconsistent group info, due to melt-operation
    pmdAnnoStatSimp = keepConsistentGroup(pmdAnnoStatSimp,'pmdGroup','lmrGroup')

    ## the same treatments apply to lmr
    pmdAnnoStatSimp = pmdAnnoStatSimp.melt(id_vars=set(pmdAnnoStatSimp.columns).difference(umr), 
                                var_name="umrGroup", value_name = 'umr')
    pmdAnnoStatSimp.loc[:,'umrGroup'] = pmdAnnoStatSimp['umrGroup'].str[0:1]
    ## drop inconsistent group info, due to melt-operation
    pmdAnnoStatSimp = keepConsistentGroup(pmdAnnoStatSimp,'pmdGroup','umrGroup')


    # Aim to exclude ambiguous annotation:  hmd-lmr or all are True, 
    # -> drop overlap regions:  isna().any
    # -> but results in insufficient genes to model
    originRows = pmdAnnoStatSimp.shape[0]
    ## 1st run  exclusion 
    #pmdAnnoStatSimp = pmdAnnoStatSimp.loc[pmdAnnoStatSimp.loc[:,['lmr','umr']].isna().any(axis=1)]
    #pmdAnnoStatSimp = pmdAnnoStatSimp.loc[pmdAnnoStatSimp.loc[:,['pmd','umr']].isna().any(axis=1)]
    #pmdAnnoStatSimp = pmdAnnoStatSimp.loc[pmdAnnoStatSimp.loc[:,['pmd','lmr']].isna().any(axis=1)]
    #print('drop any 2 overlap  reduce {:03.2f} items,before {}'.format((originRows - pmdAnnoStatSimp.shape[0])/originRows, originRows))
    ## 2nd run  exclusion 
    #pmdAnnoStatSimp = pmdAnnoStatSimp[pmdAnnoStatSimp.loc[:,['pmd','lmr','umr']].isna().any(axis=1)]
    #print('drop 3 overlap reduce {:03.2f} items,now {}'.format((originRows - pmdAnnoStatSimp.shape[0])/originRows, pmdAnnoStatSimp.shape[0]))

    # supplement the NA-NA-NA as HMD
    # supplement the T-T-T as AMD  ambiguous meathylation domain
    print(pmdAnnoStatSimp.shape)
    pmdAnnoStatSimp['hmd'] = False
    pmdAnnoStatSimp['amd'] = False
    pmdAnnoStatSimp = pmdAnnoStatSimp.reset_index()
    for index,row in enumerate(zip(pmdAnnoStatSimp['pmd'],pmdAnnoStatSimp['lmr'],pmdAnnoStatSimp['umr'])):
        #print(index)
        if pd.Series(row).isna().all():  # equal ,may faster
            pmdAnnoStatSimp.at[index, 'hmd'] = True
            pmdAnnoStatSimp.at[index, 'pmd'] = False
            pmdAnnoStatSimp.at[index, 'lmr'] = False
            pmdAnnoStatSimp.at[index, 'umr'] = False
            #print(pmdAnnoStatSimp.loc[index, :])
        elif set(row) == {True,}:    
            pmdAnnoStatSimp.at[index, 'amd'] = True  #ambiguous meathylation region
            pmdAnnoStatSimp.at[index, 'pmd'] = False
            pmdAnnoStatSimp.at[index, 'lmr'] = False
            pmdAnnoStatSimp.at[index, 'umr'] = False
            #print(pmdAnnoStatSimp.loc[index, :])
    print(pmdAnnoStatSimp.shape)
    #print(pmdAnnoStatSimp.loc[pmdAnnoStatSimp.loc[:,['amd','pmd']].all(1)].head())  amd ==TRUE & pmd ==TRUE
    
    # flatten the meanMethy
    pmdAnnoStatSimp = pmdAnnoStatSimp.melt(id_vars=set(pmdAnnoStatSimp.columns).difference(sampleName), 
                var_name="sample", value_name = 'meanMethy')
    pmdAnnoStatSimp.loc[:,'group'] = pmdAnnoStatSimp['sample'].str[0:1]
    # should match:  each group of mean methylation --- each group of methylation type 
    print(pmdAnnoStatSimp.shape[0])
    pmdAnnoStatSimp = keepConsistentGroup(pmdAnnoStatSimp,'group','pmdGroup')

    # summary types of Methy
    #flatten the methylation types
    #print(pmdAnnoStatSimp.head(30))
    pmdAnnoStatSimp = pmdAnnoStatSimp.reset_index()
    pmdAnnoStatSimp = pmdAnnoStatSimp.melt(id_vars=set(pmdAnnoStatSimp.columns).difference(('amd','hmd','pmd','lmr','umr')), 
                                var_name="typeMethy", value_name = 'type').dropna(subset=['type'])
    pmdAnnoStatSimp = pmdAnnoStatSimp.loc[pmdAnnoStatSimp.type==True,:]
    #group-by and summary
    pmdAnno_grouped = pmdAnnoStatSimp.groupby(['geneId','group','sample','cateRange','typeMethy']) #
    meanMethyDF = pmdAnno_grouped['meanMethy'].agg([np.median, np.mean, np.std])
    #print(meanMethyDF[~meanMethyDF.isna().any(axis=1)].reset_index())
    try:
        pmdAnno_grouped['meanMethy'].describe()
    except Exception as e:
        print(df_path)
        print(pmdAnno_grouped.head())
        print(e)
    else:
        #export summary table
        #pmdAnno_grouped['meanMethy'].describe().to_csv(os.path.join(out_dir,"{}meanMethydescribe_individul.txt".format(df_path)),sep='\t',index=True)
        #meanMethyDF[~meanMethyDF.isna().any(axis=1)].reset_index().to_csv(os.path.join(out_dir,"{}meanMethyDrop_individul.txt".format(df_path)),sep='\t',index=True)
        print(mem_usg(pmdAnnoStatSimp))

    # basic 没必要
    #pmdAnnoStat = pmdAnnoStat.melt(id_vars=set(pmdAnnoStat.columns).difference(pmdType), 
    #                               var_name="typeMethy", value_name = 'pmdType')
    #pmdAnnoStat = pmdAnnoStat.loc[pmdAnnoStat.pmdType == True]

    # ##diff
    # pmdAnno_groupediff = pmdAnnoStat[list(diffMethy.union(('geneId',)))]  #,'cateRange','typeMethy'
    # pmdAnno_groupediff = pmdAnno_groupediff.melt(id_vars=['geneId'], var_name="group", value_name = 'diffMethy')
    # diffMethyDF = pmdAnno_groupediff.groupby(['geneId',"group"]).agg([np.mean])
    # diffMethyDF.columns = diffMethyDF.columns.droplevel(0)
    # diffMethyDF = diffMethyDF[~diffMethyDF.isna().any(axis=1)].reset_index()
    # diffMethyDF['group'] = diffMethyDF['group'].str[0:2]
    # diffMethyDF.to_csv(os.path.join(out_dir,"{}diffMethyDrop.txt".format(df_path)),sep='\t',index=False)

    # #areaStat
    # pmdAnno_groupediff = pmdAnnoStat[list(areaStat.union(('geneId',)))]  #,'cateRange','typeMethy'
    # pmdAnno_groupediff = pmdAnno_groupediff.melt(id_vars=['geneId'], var_name="group", value_name = 'areaStat')
    # diffMethyDF = pmdAnno_groupediff.groupby(['geneId',"group"]).agg([np.mean])
    # diffMethyDF.columns = diffMethyDF.columns.droplevel(0)
    # diffMethyDF = diffMethyDF[~diffMethyDF.isna().any(axis=1)].reset_index()
    # diffMethyDF['group'] = diffMethyDF['group'].str[0:2]
    # diffMethyDF.to_csv(os.path.join(out_dir,"{}areaStatDrop.txt".format(df_path)),sep='\t',index=False)
    # print(mem_usg(pmdAnnoStat))


if __name__ == '__main__':
    main()
