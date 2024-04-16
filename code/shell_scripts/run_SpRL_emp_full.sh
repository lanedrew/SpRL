#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --time=166:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=Emp_full
#SBATCH --cpus-per-task=15


#
# modules
#
module purge
module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

#
# run scripts
#

Rscript ./code/two_stage_model/two_stage_emp_full.R

echo "All jobs successfully submitted!"