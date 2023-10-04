<!--
 * @Authour: utao.cao
-->
# epi-genetic sheep muscle development

This repository provided the workflows used in *Genome-wide epigenetic dynamics during postnatal skeletal muscle growth in Hu sheep* covering steps from reads alignment to results needed by downstream R packages.

## directory structure

- bash code has been placed into directory named by their functions.
- perl script used to filter nucleosome-free regions (NFR) reads by length before macs2 call peak
  -  bam -> bed -> filter NFR(perl) -> macs2 call peak
  -  `reads_filter.pl`
- python script used to calculate mean methylation across whole genome
  - bam -> bismark cov file -> mean methylation (python) -> quantitative matrix
  - `script/wgbs/python_WGBS/`

```bash
├── README.md
├── config
│   └── BSgenome.Oaries.NCBI.sheep-seed
├── data
│   ├── CpGislands_ovis2015.txt
│   ├── SRR_Acc_List.txt
│   ├── homerknownResults.txt
│   └── ovis2015.chrom.sizes
└── scripts
    ├── atac
    ├── motif
    ├── rnaseq
    ├── wgbs
    ├── wgs 
    ├── clipBED_autoGenerateAWKScript.sh
    ├── reads_filter.pl
    └── script_revertChromosome.awk  
```


