# scripts for downstream analysis 
After **Bsmark** outputs coverage files for all samples, we did following process:
1. pre-process:
    1. chromosome versions transitions (chr1/CMXXX/1 awk script and a mapping table )
    2. coverage calculations ( M/N, awk command ) 
2. calculate correlation between samples  **100bp as a segment**
    1. first **merge samples** by position (new detected methylated sites will recorded as 0 in non-detected samples, so if necessary, we can excluded sites by 0 numbers, **default**: not excluded)
    2. calculate sample correlation
3. Features annotation
    1. gene ID annotation by **ChIPseeker**, according to coordination;
        1. self-construct TxDb objects(from GFF3 files) 
    2. methylation type by **MethylSeekR** , utilizing HMM algorithm (**PMD, UMR, LMR**),  requiring CpG island infomation and SNP information;
    3. differential methylated sites/regions by **DSS** (ignored for next step)
4. calculate mean methylation of different types/positions for each gens
    1. dissect merged files for parallel
    2. each subfile subjected to: **pivotSummary_parallel.py**  and **pivotTableMerge.py** to collect calculated files


