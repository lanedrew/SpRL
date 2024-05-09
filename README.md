# SpRL-Code
This appendix contains the code for the paper "A Bayesian Record Linkage Approach to Tree Demography Using Overlapping LiDAR Scans".

The folders contained in this appendix are:

1. 1_simulation_study
2. 2_empirical_analysis
3. 3_figures_and_tables
4. resources

The folders and code within each folder are provided in the order they should be run. The resources folder contains supporting code (R, C++, STAN) and a placeholder folder for the empirical data used in this analysis. The empirical data can be obtained from the following ESS-DIVE repository and should be placed into the empirical_data folder within the resources folder (see the README file in that folder for additional details).

Empirical Data: https://data.ess-dive.lbl.gov/view/ess-dive-32482a38131d613-20240503T212244619

All analyses for this project were run on an HPC system using shell scripts, which an interested user should take advantage of for parallel computing.

## Dependencies

The following is a list of packages and technologies that must be installed and where they can be found.

**R version > 4.2.2**

1. ggplot2 (CRAN)
2. dplyr (CRAN)
3. tidyr (CRAN)
4. readr (CRAN)
5. data.table (CRAN)
6. purrr (CRAN)
7. patchwork (CRAN)
8. terra (CRAN)
9. sf (CRAN)
10. Rcpp (CRAN)
11. RcppArmadillo (CRAN)
12. RcppDist (CRAN)
13. xgboost (CRAN)
14. dials (CRAN)
15. parsnip (CRAN)
16. recipes (CRAN)
17. rsample (CRAN)
18. tune (CRAN)
19. workflows (CRAN)
20. yardstick (CRAN)
21. MASS (CRAN)
22. Matrix (CRAN)
23. tmvtnorm (CRAN)
24. mvtnorm (CRAN)
25. usedist (CRAN)
26. doParallel (CRAN)
27. codetools (CRAN)
28. rstan (CRAN)
29. loo (CRAN)
30. knitr (CRAN)
31. kableExtra (CRAN)
32. stringr (CRAN)
33. ggpubr (CRAN)
34. latex2exp (CRAN)
35. gtools (CRAN)
36. mcclust (CRAN)
