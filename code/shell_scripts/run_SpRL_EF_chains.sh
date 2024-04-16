#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=Launch_Emp_Full

for chain in 1 2 3 4
do
for sigma2_prior in "weak"
do

sbatch --job-name=EF_Chain_${chain}_${sigma2_prior} --account=csu68_alpine1 --qos=long --partition=amilan --time=167:50:00 --mem=64G --export=chain=$chain,sigma2_prior=$sigma2_prior ./code/shell_scripts/call_SpRL_emp_full.sh

done
done

echo "All jobs successfully submitted!"