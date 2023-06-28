#!/bin/bash
#SBATCH --time=1000:00:00  # walltime limit (HH:MM:SS)
#SBATCH --cpus-per-task=20 #number of threads per task
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"

# HISAT2 index of oyster genome
####################################################################

echo "******************** start script" $(date)

module load HISAT2/2.2.1-gompi-2021b
module list

# Setup
####################################################################

threads=16
# genome downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/
INPUT=/data/marine_diseases_lab/shared/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_genomic.fna
INDEX=/data/marine_diseases_lab/shared/cvir_hisat2.2.1index/cvirhisat

# Analysis
####################################################################

echo "******************** processing" $INPUT

hisat2-build -p $threads -f $INPUT $INDEX

echo "******************** done script" $(date)
