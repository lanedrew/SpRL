#!/bin/bash

module purge
module load intel/2022.1.2
module load R/4.2.2
module load gdal/3.5.0 geos/3.10.2 proj/8.2.1 sqlite/3.38.01 cmake/3.25.0

case "$model" in
    #case 1
    "pe") Rscript ./code/two_stage_model/two_stage_growth_pe_cpp.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL} ${n_thresh:-NULL};;

    #case 2
    "la") Rscript ./code/two_stage_model/two_stage_growth_la_cpp.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL} ${n_thresh:-NULL} ${link_set:-NULL};;
    
    #case 3
    "ndm") Rscript ./code/naive_growth_models/ndm_growth_sim.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL};;
    
    #case 4
    "TL") Rscript ./code/two_stage_model/growth_sim_results_TL.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL};;
    
    #case 5
    "RL") Rscript ./code/two_stage_model/two_stage_sim_cpp.R ${density:-NULL} ${noise:-NULL} ${alpha:-NULL} ${index:-NULL} ${n_thresh:-NULL};; 
esac

