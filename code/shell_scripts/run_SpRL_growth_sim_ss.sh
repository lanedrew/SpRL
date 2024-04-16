#!/bin/bash

#SBATCH --account=csu68_alpine1
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=14:00:00
#SBATCH --partition=amilan

ncores=6

function check_active_jobs {
    num_active_jobs=$(squeue -u $USER | grep -c "R\|PD")
    echo "$num_active_jobs"
}

MAX_TOTAL_JOBS=999    # Maximum total number of active jobs


#
# modules
#
module purge

#
# run scripts
#

for model in "pe"
do
for density in "low" "med" "high"
do
for noise in "small" "medium" "large"
do
for alpha in 1 2 3
do
for index in $(seq 1 100)
do

while [ $(check_active_jobs) -ge $MAX_TOTAL_JOBS ]; do
    sleep 60  # Sleep for a minute before checking again
done

sbatch --account=csu68_alpine1 --job-name=twost_${model}_${density}_${noise}_${alpha}_${index} --qos=normal --partition=amilan --time=14:00:00 --cpus-per-task=$ncores --export=density=$density,noise=$noise,alpha=$alpha,index=$index,model=$model ./code/shell_scripts/call_SpRL_growth_two_stage.sh
done
done
done
done
done


# for density in "low" "med" "high"
# do
# for noise in "medium"
# do
# for alpha in 2
# do
# for index in 1 2 3
# do
# 
# while [ $(check_active_jobs) -ge $MAX_TOTAL_JOBS ]; do
#     sleep 60  # Sleep for a minute before checking again
# done
# 
# sbatch --job-name=MS_${density}_${noise}_${index} --qos=normal --partition=amilan --time=24:00:00 --cpus-per-task=$ncores --export=density=$density,noise=$noise,alpha=$alpha,index=$index ./code/shell_scripts/call_SpRL_growth_sim_study.sh
# done
# done
# done
# done


# for density in "med" "high"
# do
# for noise in "small" "medium" "large"
# do
# for alpha in 1 2 3
# do
# for index in $(seq 1 100)
# do
# 
# while [ $(check_active_jobs) -ge $MAX_TOTAL_JOBS ]; do
#     sleep 60  # Sleep for a minute before checking again
# done
# 
# sbatch --account=csu68_alpine1 --job-name=twost_${density}_${noise}_${alpha}_${index} --qos=normal --partition=amilan --time=24:00:00 --cpus-per-task=$ncores --export=density=$density,noise=$noise,alpha=$alpha,index=$index ./code/shell_scripts/call_SpRL_growth_two_stage.sh
# done
# done
# done
# done



# for density in "high"
# do
# for noise in "small" "medium" "large"
# do
# for alpha in 1 2 3
# do
# for index in $(seq 1 100)
# do
# 
# # while [ $(check_active_jobs) -ge $MAX_TOTAL_JOBS ]; do
# #     sleep 60  # Sleep for a minute before checking again
# # done
# 
# sbatch --job-name=jm_${density}_${noise}_${alpha}_${index} --qos=normal --partition=amilan --time=24:00:00 --cpus-per-task=$ncores --export=density=$density,noise=$noise,alpha=$alpha,index=$index ./code/shell_scripts/call_SpRL_growth_sim_study.sh
# done
# done
# echo "Noise level: $noise submitted."
# done
# done



echo "All jobs successfully submitted!"