## 16SampliconAnalysis
This folder contains scripts used to perform QIIME2 analysis on the 16S rRNA amplicon data. The raw sequences generated for this study can be found in the NCBI SRA under BioProject no. PRJNA599137. See [`PRJNA599137_NCBI_16Smetatrans.xlsx`](/PRJNA599137_NCBI_16Smetatrans.xlsx) for accession numbers for each sample.  

## [scripts/](scripts/)
QIIME2 bash script and Rmd files  
- qiime2_PJallsamples.sh - main QIIME2 script with all commands
- qiime2_trainclassifier.sh - QIIME2 to train taxonomy classifier with SILVA 138-99
- PJ_16Samplicon_analysis.Rmd and knitted outputs

## [qiime2output/](qiime2output/)
QIIME2 artifacts and processed files  
- all qza and qzv files generated during the QIIME2 pipeline

## [figures/](figures/)
Output from Rmd file in scripts/  

## [metadata/](metadata/)
- **sample-manifest_PJ_V6.csv** - file to import raw sequence files into QIIME2 .qza format using `qiime tools import`
- **PJ_V6Samples_Metadata.txt** - metadata file formatted for QIIME2 import during `qiime feature-table summarize`  
