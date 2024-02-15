# nepal_earthquake
Replication file for Matthew Gordon's job market paper Visibility and Vulnerability: Aid allocations after the 2015 Nepal Earthquake. Paper available from: https://sites.google.com/view/mdgordon/.

The starting point for the analysis is the World Bank Household Risk and Vulnerability Survey, available from: https://microdata.worldbank.org/index.php/catalog/3705/related-materials.

First run the WB_HRVS_processing.R script to create full_panel.csv.

WB_HRVS_descriptive.Rmd generates the descriptive statistics, summary stats, and Townsend regressions. eq_survey_analysis.Rmd generates the correlates of earthquake damages discussed in the paper, with tables in the appendix.

WB_HRVS_estimation_pooled.Rmd generates the regression discontinuity results. rdhelpers.R and rdrobustr.r contain functions called in this script. eq_survey_analysis.r generates the tables in the appendix showing the determinants of earthquake damages.

gmmprep.r adjusts the raw data for lifecycle and household composition to prepare it for the structural estimation step. gmm.r uses GMM to calibrate the structural model. VFIfunctions.R contains helper functions used by these scripts.

VFI.Rmd uses the final calibrated value function to perform counterfactuals and compare model simulations to the RD results.
