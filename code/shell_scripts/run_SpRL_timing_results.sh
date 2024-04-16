#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=gen_timing_res

# module purge
# module load intel/2022.1.2
# module load R/4.2.2
# module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0
# 
# Rscript ./code/empirical_data/timing_results/generate_timing_results.R 

for box_margin in  "150"
do

sbatch --job-name=timing_300_${box_margin} --account=csu68_alpine1 --qos=normal --partition=amilan --time=24:00:00 --mem=16G --export=box_margin=$box_margin ./code/shell_scripts/call_SpRL_timing_results.sh

done