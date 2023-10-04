#!/usr/bin/env python
# coding: utf-8
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -J pmd
#SBATCH -e %J_err.log
#SBATCH -o %J_out.log
from multiprocessing import Pool
import sys
import os
import pandas as pd
import subprocess
#https://www.cnblogs.com/shiyanhe/p/13604707.html  Python创建大量线程时遇上OpenBLAS blas_thread_init报错怎么办？
os.environ['OPENBLAS_NUM_THREADS'] = '2'

def job1(k):
    print('calling >>> python WGBS_PMD_reanalysis_chr.py %s %d  <<<' % (k[0], k[1]))
    #subprocess.run(' python -u WGBS_PMD_reanalysis_chr.py %s %d  500 > log_500bp 2>&1 ' % (k[0], k[1]),shell=True )
    #subprocess.run(' python -u WGBS_PMD_reanalysis_chr.py %s %d  100 > log_100bp 2>&1 ' % (k[0], k[1]),shell=True  )
    subprocess.run(' python -u WGBS_PMD_reanalysis_chr.py %s %d  300 > log_300bp 2>&1 ' % (k[0], k[1]),shell=True  )
    subprocess.run(' python -u WGBS_PMD_reanalysis_chr.py %s %d  1000 > log_1000bp 2>&1 ' % (k[0], k[1]),shell=True  )
    subprocess.run('free -h',shell=True)

def job2(k):
    print('calling >>> python remodify.py %s %d  1000 <<<' % (k[0], k[1]))
    #subprocess.run('python -u remodify.py %s %d 300' % (k[0], k[1]),shell=True  )
    subprocess.run('python -u remodify.py %s %d 100' % (k[0], k[1]),shell=True  )
    subprocess.run('free -h',shell=True)
    #subprocess.run('python -u remodify.py %s %d 1000' % (k[0], k[1]),shell=True  )

autoChromosome = pd.read_csv("GCA_ovis4.02015_report.txt",
                             sep="\t",nrows=26)

chromo = {}
for row in zip(autoChromosome["Sequence-Length"],autoChromosome["UCSC-style-name"]):
    chromo[row[1]] = row[0]

pool = Pool()
data_list =  [i for i in chromo.items() ]
pool.map(job2, data_list)
#pool.map(job1, data_list)
pool.close()
pool.join()
