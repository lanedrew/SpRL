#!/bin/bash
ncores=1

#SBATCH --job-name=SpRL_data_gen
#SBATCH --qos=normal
#SBATCH --partition=smem
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=$ncores
#SBATCH --mem=20G

module purge
source /curc/sw/anaconda3/latest
conda activate Renv

Rscript code/datasets_growth_sim_study.R
