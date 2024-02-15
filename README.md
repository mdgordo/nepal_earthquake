# nepal_earthquake
Replication file for Matthew Gordon's job market paper Visibility and Vulnerability: Aid allocations after the 2015 Nepal Earthquake. Paper available from: https://sites.google.com/view/mdgordon/.

The starting point for the analysis is the World Bank Household Risk and Vulnerability Survey, available from: https://microdata.worldbank.org/index.php/catalog/3705/related-materials.

First run the WB_HRVS_processing.R script to create full_panel.csv.

WB_HRVS_descriptive.Rmd generates the descriptive statistics, summary stats, and Townsend regressions. eq_survey_analysis.Rmd generates the correlates of earthquake damages discussed in the paper, with tables in the appendix.

WB_HRVS_estimation_pooled.Rmd generates the regression discontinuity results. rdhelpers.R contains functions called in this script.

gmmprep.R prepares the raw data for the structural calibration. gmm.R calibrates the model. VFIfunctions.R contains helper functions used by both of these scripts. 

VFI.Rmd uses the final calibrated value function to perform counterfactuals and compare model simulations to the RD results.
