
<!-- README.md is generated from README.Rmd. Please edit that file -->
S2MET
=====

Description
-----------

Project repository for the spring two-row (**S2**) barley multi-environment trial (**MET**) project. This project will look at the impact of genotype-by-environment interaction when performing genomewide prediction. The project includes a panel of 183 training population lines, 50 validation population lines, and 9 repeated checks.

Script Order
------------

This is the order in which scripts should be executed to replicate the results. All scripts are located in the [Scripts](https://github.com/neyhartj/S2MET/tree/master/Scripts) directory.

Phenotypic data
1. `phenotype_data_adjustment.R` - Calculate genotypic means and calculate heritability in each environment.
2. `phenotype_data_summary.R` - Summarize phenotypic data across all environments, including distributions and heritability.

Environmental variables
1. `collect_environmental_variables.R` - Collect environmental covariate data from NOAA and Soil Survey databases.
2. `manipulate_analyze_environmental_variables.R` - Calculate summary variables, including growing degree days (GDD), accumulated GDD (AGDD), and photothermality.
3. `model_environmental_variables.R` - Find the variables that are highly correlated with the environmental mean and create covariance and distance matrices.

Clustering
1. `environmental_clustering.R` - Calculate great circle distance and phenotypic distance, then create hierarchical clusters based on the distance objects.
2. `cluster_heritability.R` -
