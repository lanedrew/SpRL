#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=1:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=scrps_model_comp
#SBATCH --mem=64GB

module purge
module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

# Rscript ./code/results_processing/all_model_comp_scrps.R
Rscript ./code/results_processing/all_model_comp_scrps_LA.R