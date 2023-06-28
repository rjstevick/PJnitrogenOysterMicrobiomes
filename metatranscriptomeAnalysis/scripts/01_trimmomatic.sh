#!/bin/bash

#SBATCH --time=1000:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=20 #number of threads per task
#SBATCH --array=1-36%18 #number of total jobs(based on number of files)%number of jobs at a time
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"

# Trimmomatic quality control
####################################################################

echo "******************** start script" $(date)
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load Trimmomatic
module list


# Setup
####################################################################

threads=16
trimlog=/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/01_QC/trimlog

INPUT_DIR=/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/00_RAW
OUTPUT_DIR=/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/01_QC

# get file names
f1=`ls $INPUT_DIR/*R1_001.fastq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
base=`echo $f1 | sed 's/_L001_R1_001.fastq.gz//g' | sed 's+/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/00_RAW/++g'`
f2=`echo $f1 | awk -F "R1" '{print $1 "R2" $2}'`


# Analysis
####################################################################

echo "******************** processing" $base

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $threads -phred33 \
    -trimlog $trimlog.$base".txt" \
    $f1 $f2 \
    $OUTPUT_DIR/$base".R1.clean.fastq.gz" $OUTPUT_DIR/$base".R1_unpaired.clean.fastq.gz" \
    $OUTPUT_DIR/$base".R2.clean.fastq.gz" $OUTPUT_DIR/$base".R2_unpaired.clean.fastq.gz" \
    ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15:keepBothReads MINLEN:140 \

echo "******************** done script" $(date)
