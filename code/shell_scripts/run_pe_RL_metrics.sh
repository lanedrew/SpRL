#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=6:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=pe_RL_metrics


for density in "low" "med" "high"
do
for noise in "small" "medium" "large"
do
for alpha in 1 2 3
do

sbatch --account=csu68_alpine1 --job-name=pe_RL_mets_${density}_${noise}_${alpha} --qos=normal --partition=amilan --time=7:00:00 --cpus-per-task=15 --export=density=$density,noise=$noise,alpha=$alpha ./code/shell_scripts/call_pe_RL_metrics.sh

done
done
done

echo "All jobs successfully submitted!"