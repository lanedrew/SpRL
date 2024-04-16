#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=Launch_Emp_Full

for chain in 3
do

sbatch --job-name=EF_Chain_ext${chain} --account=csu68_alpine1 --qos=long --partition=amilan --time=167:50:00 --cpus-per-task=15 --export=chain=$chain ./code/shell_scripts/call_SpRL_EF_ext.sh

done

echo "All jobs successfully submitted!"