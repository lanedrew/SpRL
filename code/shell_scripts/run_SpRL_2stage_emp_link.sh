#!/bin/sh
ncores=3

#
# modules
#
module purge

#
# run scripts
#
for subset in "1" "2" "3" "4"
do
sbatch --job-name=two_stage_full --qos=long --partition=amilan --time=166:00:00 --cpus-per-task=$ncores --export=subset=$subset call_SpRL_2stage_emp_link.sh
done
