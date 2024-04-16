#!/bin/bash

#SBATCH --account=csu68_alpine1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --job-name=emp_pooling
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=24:00:00

module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

Rscript ./code/results_processing/empirical_chains_pooling.R 
