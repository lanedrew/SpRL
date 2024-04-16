#!/bin/bash

module purge
module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

Rscript ./code/two_stage_model/two_stage_emp_full_la.R ${covars:-NULL} ${index:-NULL} ${mort_threshold:-NULL} ${model_type:-NULL}
# Rscript ./code/results_processing/calculate_scrps_LA.R ${covars:-NULL} ${index:-NULL} ${mort_threshold:-NULL} ${model_type:-NULL}
# Rscript ./code/naive_growth_models/polygon_overlap_matching_growth_model.R ${covars:-NULL} ${index:-NULL} ${mort_threshold:-NULL}
# Rscript ./code/naive_growth_models/near_distance_matching_growth_model.R ${covars:-NULL} ${index:-NULL} ${mort_threshold:-NULL}