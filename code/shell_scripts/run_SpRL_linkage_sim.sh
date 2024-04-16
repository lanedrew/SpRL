#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=Launch_RL_sim_study


function check_active_jobs {
    num_active_jobs=$(squeue -u $USER | grep -c "R\|PD")
    echo "$num_active_jobs"
}

MAX_TOTAL_JOBS=999    # Maximum total number of active jobs


for model in "RL"
do
for density in "low" "med" "high"
do
for noise in "small" "medium" "large"
do
for alpha in 1 2 3
do
for index in $(seq 1 100)
do
for n_thresh in 25
do

while [ $(check_active_jobs) -ge $MAX_TOTAL_JOBS ]; do
    sleep 60  # Sleep for a minute before checking again
done

sbatch --job-name=twost_${model}_${density}_${noise}_${alpha}_${index} --qos=normal --partition=amilan --time=1:00:00 --cpus-per-task=4 --export=density=$density,noise=$noise,alpha=$alpha,index=$index,model=$model,n_thresh=$n_thresh ./code/shell_scripts/call_SpRL_growth_two_stage.sh
done
done
done
done
done
done

