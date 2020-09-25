#!/bin/bash

# #PBS -l walltime=48:00:00,mem=62gb,nodes=1:ppn=8
#PBS -l walltime=08:00:00,mem=62gb,nodes=1:ppn=12
# #PBS -N loeo_predictions_fr
#PBS -N lolo_predictions_fr
# #PBS -N cv_7525_predictions_fr
# #PBS -N lolo_predictions_longterm_fr
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/Scripts/Predictions/

module load R/3.5.2_mkl

## Factorial regression 

# Leave-one-environment-out
# Rscript environment_loeo_predictions2.R

# # Leave-one-location-out
Rscript environment_lolo_predictions2.R


## 75-25 cross-validation
# Rscript cross_validation_7525_predictions2.R


## Locations predictions using long-term covariates
# Rscript environment_lolo_longterm_covariates_predictions2.R
