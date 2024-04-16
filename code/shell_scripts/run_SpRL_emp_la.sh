#!/bin/bash

#SBATCH --account=csu68_alpine1 
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=1:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=Launch_Emp_LA_Full


function check_active_jobs {
    num_active_jobs=$(squeue -u $USER | grep -c "R\|PD")
    echo "$num_active_jobs"
}

MAX_TOTAL_JOBS=199    # Maximum total number of active jobs


for covars in "all"
do
for index in $(seq 1 100)
do
for mort_threshold in 90
do
for model_type in "skew_t" "skew_normal" "normal" "MLR"
do

# while [ $(check_active_jobs) -ge $MAX_TOTAL_JOBS ]; do
#     sleep 60  # Sleep for a minute before checking again
# done

# sbatch --job-name=emp_LA_${covars}_${index}_${mort_threshold} --account=csu68_alpine1 --qos=normal --partition=amilan --time=24:00:00 --cpus-per-task=6 --export=covars=$covars,index=$index,mort_threshold=$mort_threshold ./code/shell_scripts/call_SpRL_emp_la.sh
sbatch --job-name=emp_LA_${model_type}_${index} --account=csu68_alpine1 --qos=normal --partition=amilan --time=24:00:00 --mem=32G --export=covars=$covars,index=$index,mort_threshold=$mort_threshold,model_type=$model_type ./code/shell_scripts/call_SpRL_emp_la.sh
# sbatch --job-name=ndm_${covars}_${mort_threshold} --account=csu68_alpine1 --qos=normal --partition=amilan --time=24:00:00 --mem=32G --export=covars=$covars,index=$index,mort_threshold=$mort_threshold ./code/shell_scripts/call_SpRL_emp_la.sh

done
done
done
done

echo "All jobs successfully submitted!"