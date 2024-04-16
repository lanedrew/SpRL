#!/bin/bash

source /curc/sw/anaconda3/latest
conda activate Renv_new

Rscript ./code/two_stage_model/two_stage_sim_cpp.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL}