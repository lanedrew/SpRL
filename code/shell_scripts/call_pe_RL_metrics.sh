#!/bin/bash

module purge
module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

# Rscript ./code/results_processing/two_stage_pe_RL_metrics.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL}
Rscript ./code/results_processing/two_stage_la_RL_metrics.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL}
# Rscript ./code/results_processing/ndm_RL_metrics.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL}
