#!/bin/bash
#SBATCH --time=1000:00:00  # walltime limit (HH:MM:SS)
#SBATCH --cpus-per-task=16 #number of threads per task
#SBATCH --array=1-36%18 #number of total jobs(based on number of files)%number of jobs at a time
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"

# HISAT2 mapping to oyster genome
####################################################################

echo "******************** start script" $(date)
echo "...SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "...SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load HISAT2/2.2.1-gompi-2021b
module list


# Setup
####################################################################

threads=16

INPUT_DIR=/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/01_QC
OUTPUT_DIR=/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/02_MAPPING_oyster
INDEX=/data/marine_diseases_lab/shared/cvir_hisat2.2.1index/cvirhisat

# get file names
f1=`ls $INPUT_DIR/*R1.clean.fastq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
base=`echo $f1 | sed 's+.R1.clean.fastq.gz++g' | sed 's+/data/marine_diseases_lab/rebecca/PJsequencing/metatranscriptomes/01_QC/++g'`
f2=`echo $f1 | awk -F "R1" '{print $1 "R2" $2}'`


# Analysis
####################################################################

echo "******************** processing" $base

hisat2 -t \ # print wall-clock time taken by search phases
    -p $threads \ # number of alignment threads to launch
    -x $INDEX \ # Index filename prefix (minus trailing .X.ht2)
    -1 $f1 -2 $f2 \ # Files with mates, paired with files
    -S $base.sam \ # File for SAM output
    --un-conc $base.pairedunmap.fastq \ # write pairs that didn't align concordantly to <path>
    --met-file $base.log # send metrics to file at <path>

# move all files to the output directory
mv $base.* $OUTPUT_DIR

echo "******************** done script" $(date)
