
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BarleyEnvironmentGWS - Scripts

Below is a description of scripts used in this analysis. Within each
subfolder, script names have a numeric prefix listing the order in which
the scripts should be run. The order of the subfolders below also
corresponds to the order of script execution.

## PhenotypicAnalysis

1.  `01_phenotypic_data_summary.R` - Summarizes phenotypic data and
    estimates variance components under different models.

## EnvironmentalCovariables

1.  `01_query_environmental_data.R` - retrieves raw daily weather data
    and soil data from databases.
2.  `02_run_crop_growth_model.R` - uses the daily weather data to
    parameterize crop growth models.
3.  `03_summarize_environmenta_data_by_growth_stage.R` - summarizes the
    raw weather data according to the growth stages predicted by the
    crop growth model.
4.  `04_define_environmental_covariables.R` - defines environmental
    covariates and runs diagnostics on the covariates.
5.  `05_covariate_phenotypic_modeling.R` - performs feature selection or
    assigns feature importance to the covariates.
6.  `06_covariate_phenotypic_modeling_location.R` - same as \#5, but
    using the historical weather data.
7.  `07_interval_covariate_phenotypic_modeling.R`- same as \#5, but
    using environmental covariates defined by arbitrary 20-day windows
    instead of predicted growth stages.
8.  `08_interval_covariate_phenotypic_modeling_location.R` - same as
    \#6, but using environmental covariates defined by arbitrary 20-day
    windows instead of predicted growth stages.

## Predictions

1.  `01_environmental_predictions.R` - runs leave-one-environment-out
    cross-validation using the concurrent weather and soil data.
2.  `02_location_predictions.R` - runs leave-one-location-out
    cross-validation using historical weather and soil data.
3.  `03_prediction_analysis.R` - analyzes predictions, calculates
    predictive ability and prediction accuracy.
4.  `A01_environmental_predictions_varcomp.R` - same as \#1, but using
    REML estimates of variance components from the
    `01_phenotypic_data_summary.R` script.
5.  `A02_location_predictions_varcomp.R` - same as \#2, but using REML
    estimates of variance components from the
    `01_phenotypic_data_summary.R` script.
6.  `A03_varcomp_prediction_analysis.R` - same as \#3, but using REML
    estimates of variance components from the
    `01_phenotypic_data_summary.R` script.
