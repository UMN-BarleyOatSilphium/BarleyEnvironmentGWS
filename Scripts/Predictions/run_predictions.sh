#!/bin/bash

#PBS -l walltime=36:00:00,mem=62gb,nodes=1:ppn=12
# #PBS -l walltime=02:00:00,mem=62gb,nodes=1:ppn=12
# #PBS -N loeo_predictions_fr
#PBS -N lolo_predictions_fr
# #PBS -N loo_predictions_fr
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

# # Leave-one-out from individual covariate selection
# Rscript environment_loo_predictions2_fr_samples.R


# # Leave-one-out and external validation
# Rscript environment_loo_predictions2.R


