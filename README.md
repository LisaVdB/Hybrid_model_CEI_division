# Hybrid_model_CEI_division
Exchange of molecular and cellular information: a hybrid model that integrates stem cell divisions and key regulatory interactions.

The following three files were used to identify sensitive parameters within the hybrid model:
- gif1_sensitivity_analysis.m
- gif1_dy.m (with the ODEs used in the hybrid model)
- run_SA.m (to run the sensitivity analysis, returning S_results and ST_results)
The ST_results were used to calculate the global sensitivity index.

The following two files allowed for the estimation of 11 parameters gif1_param_est.m file:
- gif1_param_est.m
- gif1_dy_est.m (with the ODEs used in the hybrid model)

A second hybrid model was constructed to accurately fit an3 mutant data. For this 6 parameters were re-estimated with the following files:
- gif1_param_est_mutant.m
- gif1_dy_est_mutant.m (with the ODEs used in the hybrid model)

The hybrid models were created in SimBiology to model, simulate, and analyze the dynamic systems. The following two .sbproj files contain the first and second model:
- GIF1_Final.sbproj
- GIF1_an3_Final.sbproj (an additional factor X and its upstream and downstream regulations were added)
