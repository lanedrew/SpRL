#!/bin/bash

#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --job-name=alpine_plots
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=6:00:00

module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

Rscript ./code/results_processing/la_plots_alpine.R 
