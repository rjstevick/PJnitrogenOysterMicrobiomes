# Metatranscriptome analysis

## Goals

- ??
- Analyze expression of bacterial nitrogen cycling genes upon nutrient enrichment in oyster gut, inner shell, and outer shell swab samples.

## In this README

- [Contents](#contents)
- [Bioinformatic Analysis](#bioinformatic-analysis)
- [Results Summary](#results-summary)


---------------------------------------
<br/>

# Contents

## [output/](output/)
This folder contains processed files and summaries to generate plots.  

## [figures/](figures/)
This folder contains all generated figure files from [scripts/](scripts/)

## [scripts/](scripts/)
- All bash scripts used to analyze data with slurm scheduling framework

## [metadata/](metadata/)
- Tables of metadata from environmental monitoring

---------------------------------------
<br/>

# Bioinformatic Analysis


## 0. Check Data

**FastQC**

```bash
fastqc $INPUT_DIR/*gz -o $OUTPUT_DIR
multiqc $OUTPUT_DIR/*
```


## 1. Quality Control

**Trimmomatic** - scripts/01_trimmomatic.sh  

```bash
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $threads -phred33 \
    -trimlog $trimlog.$base".txt" \
    $f1 $f2 \
    $OUTPUT_DIR/$base".R1.clean.fastq.gz" $OUTPUT_DIR/$base".R1_unpaired.clean.fastq.gz" \
    $OUTPUT_DIR/$base".R2.clean.fastq.gz" $OUTPUT_DIR/$base".R2_unpaired.clean.fastq.gz" \
    ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15:keepBothReads MINLEN:140 \
```

**Output:** 01_QC/*.clean.fastq.gz


## 2. Mapping to oyster genome

Use assembly 3.0 GCF_002022765.2_C_virginica-3.0_genomic.fna  
Index made with scripts/02_indexhisat2.sh

**HISAT2 v2.2.1** - scripts/02_hisat2.sh  

```bash
hisat2 -t \ # print wall-clock time taken by search phases
    -p $threads \ # number of alignment threads to launch
    -x $INDEX \ # Index filename prefix (minus trailing .X.ht2)
    -1 $f1 -2 $f2 \ # Files with mates, paired with files
    -S $base.sam \ # File for SAM output
    --un-conc $base.pairedunmap.fastq \ # write pairs that didn't align concordantly to <path>
    --met-file $base.log # send metrics to file at <path>
```

**Output:** 02_MAPPING_oyster/*.sam
