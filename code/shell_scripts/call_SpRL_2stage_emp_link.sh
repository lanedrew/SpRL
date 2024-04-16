#!/bin/bash

source /curc/sw/anaconda3/latest
conda activate Lane-env

Rscript ./code/SpRL_2stage_linkage_full.R ${subset:-NULL}
