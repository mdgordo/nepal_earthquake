# nepal_earthquake
Replication file for Targeting Disaster Aid: A Structural Evaluation of a Large Earthquake Reconstruction Program. Paper available from: https://sites.google.com/view/mdgordon/.

The starting point for the analysis is the World Bank Household Risk and Vulnerability Survey, available from: https://microdata.worldbank.org/index.php/catalog/3705/related-materials.

First run the WB_HRVS_processing.R script to create full_panel.csv.

WB_HRVS_descriptive.Rmd generates the descriptive statistics, summary stats, and Townsend regressions. eq_survey_analysis.Rmd generates the correlates of earthquake damages discussed in the paper, with tables in the appendix.

WB_HRVS_estimation.Rmd generates the regression discontinuity results. rdhelpers.R and rdrobustr.r contain functions called in this script.

gmmprep.r adjusts the raw data for lifecycle and household composition to prepare it for the structural estimation step. gmm.r uses GMM to calibrate the structural model. VFIfunctions.R contains helper functions used by these scripts.

VFI.Rmd uses the final calibrated value function to perform counterfactuals and compare model simulations to the RD results.
