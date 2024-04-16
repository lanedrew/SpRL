#!/bin/bash

# source /curc/sw/anaconda3/latest
# conda activate Renv_new

module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

# Rscript ./code/joint_model/joint_sim_cpp.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL}
Rscript ./code/joint_model/joint_sim_cpp_ms.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL}
# Rscript ./code/joint_model/joint_sim_fixed_alpha_cpp.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL}
