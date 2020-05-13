#!/bin/bash

#PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=12
# #PBS -N loeo_predictions
# #PBS -N loeo_predictions_fr
# #PBS -N loyo_predictions
# #PBS -N lolo_predictions
#PBS -N lolo_predictions_fr
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/Scripts/Predictions/

module load R/3.5.2_mkl

# # Leave-one-environment-out
# Rscript environment_loeo_predictions.R

# # Leave-one-environment-out - factorial regression
# Rscript environment_loeo_predictions2.R

# # Leave-one-year-out
# Rscript environment_loyo_predictions.R

# Leave-one-location-out
# Rscript environment_lolo_predictions.R

# Leave-one-location-out - factorial regression
Rscript environment_lolo_predictions2.R


