# Empirical Analysis
The numbered scripts in this folder can be run to recreate the empirical data analysis performed for the paper. Scripts should be run in increasing numerical order.

We note that many of these scripts are designed to be run on an HPC utilizing parallelization to maximize efficiency. For this analysis we allowed the following inputs to vary with levels:

### (Record Linkage Model)
- N Threshold: "10" "25" 
- $\sigma^2$ prior: "weak" "strong"
- Chain Index: "1" "2" "3" "4"

### (All Models)
- Mortality Threshold: "80" "90" "100"
- Model Type: "skew_t" "skew_normal" "normal" "MLR"

### (Timing Results)
- Box Margin: "1" "2" "3" "5" "10" "25" "50" "100" "150" "200" "300"

## Script Descriptions
We provide descriptions of each script in the sequence below.

1. 1_two_stage_RL_full.R: 

    This script runs the initial sampling period for the full empirical dataset given a set of input arguments (see above).
  
2. 2_two_stage_RL_full_ext.R: 

    This script runs the extended sampling period for the full empirical dataset given a set of input arguments (see above).
  
3. 3_generate_competition_metrics.R: 

    This script generates the competition metric covariates for the 2015 and 2019 datasets.
  
4. 4_pool_linkage_chains.R: 

    This script pools the results from the multiple linkage chains generated in steps 1 and 2.
  
5. 5_generate_pooled_sample_index.R: 

    This script generates the random sample index for fitting the downstream model using the LA approach.
  
6. 6_two_stage_growth_full_LA.R: 

    This script runs the downstream growth model on the empirical data given a posterior linkage sample and a set of input arguments (see above).

7. 7_generate_nearest_distance_matching_linkage.R: 

    This script generates the nearest distance matching linkage for the empirical data.
  
8. 8_nearest_distance_matching_growth_model.R

    This script runs the downstream growth model on the empirical data given the NDM linkage and a set of input arguments (see above).

9. 9_polygon_overlap_matching_growth_model.R

    This script runs the downstream growth model on the empirical data given the POM linkage and a set of input arguments (see above).

10. 10_generate_timing_results.R

    This script generates the timing results for the record linkage model fit to the empirical data given a set of input arguments (see above).