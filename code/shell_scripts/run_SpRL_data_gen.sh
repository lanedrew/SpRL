#!/bin/bash

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --partition=amilan

ncores=8

function check_active_jobs {
    num_active_jobs=$(squeue -u $USER | grep -c "R\|PD")
    echo "$num_active_jobs"
}

MAX_TOTAL_JOBS=999     # Maximum total number of active jobs


#
# modules
#
module purge


for density in "low" "med" "high"
do
for tau2 in .5
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

sbatch --job-name=dg_${density}_${alpha}_${index} --qos=normal --partition=amilan --time=24:00:00 --cpus-per-task=$ncores \
--export=density=${density},tau2=${tau2},alpha=${alpha},noise=${noise},index=${index} ./code/call_SpRL_data_gen.sh

done
done
done
done
done