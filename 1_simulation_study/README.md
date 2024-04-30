# Simulation Study Details

The numbered scripts in this folder can be run to recreate the simulation study performed for the paper. Scripts should be run in increasing numerical order.

We note that many of these scripts are designed to be run on an HPC utilizing parallelization to maximize efficiency. For this simulation study we varied the following inputs and their levels:

(All Simulations) \\
Density: "low" "med" "high"\\
Noise Level: "small" "medium" "large"\\
Alpha Value: "1" "2" "3"\\
Index: 1,...,100\\

(Linkage Averaging Additions) \\
N Threshold: "10" "25" \\
Linkage Set: "1" "2" \\

## Script Descriptions
We provide descriptions of each script in the sequence below.

1. 1_size_pred_model_build.R:\\
  This script builds the GBM models used in the data simulation procedure.
  
2. 2_emp_data_subset_processing.R:\\
  This script generates the empirical subsets used as our references for the low, medium, and high density simulation settings.
  
3. 3_dataset_generator.R:\\
  This script generates a simulated dataset given a set of input arguments (see above).
  
4. 4_two_stage_RL.R:\\
  This script runs and processes the linkage performance on a simulated dataset given a set of input arguments (see above).
  
5. 5_two_stage_growth_LA.R:\\
  This scripts runs the downstream growth model on a simulated dataset given a posterior linkage sample and a set of input arguments (see above).
  
6. 6_two_stage_growth_NDM.R:\\
  This script generates the Nearest Distance Matching linkage and runs the downstream growth model on a simulated dataset given a set of input arguments (see above).
  
7. 7_two_stage_growth_TL.R:\\
  This script runs the downstream growth model on a simulated dataset given the true linkage structure and a set of input arguments (see above).
  
