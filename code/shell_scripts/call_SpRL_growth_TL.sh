#!/bin/bash

module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

Rscript ./code/two_stage_model/growth_sim_results_TL.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL}
