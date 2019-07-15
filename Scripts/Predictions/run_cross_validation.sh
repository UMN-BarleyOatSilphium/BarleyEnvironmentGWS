#!/bin/bash

# #PBS -l walltime=48:00:00,mem=62gb,nodes=1:ppn=8
# #PBS -l walltime=72:00:00,mem=62gb,nodes=1:ppn=16
#PBS -l walltime=24:00:00,mem=64gb,nodes=1:ppn=8
# #PBS -N cross-validation-cv1
# #PBS -N cross-validation-cv2
# #PBS -N cross-validation-cv00
# #PBS -N parent-offspring-validation
#PBS -N prediction-validation-sample
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/Scripts/Predictions/

module load R/3.5.0
# module load R/3.5.2_mkl

# Cross validation
# Rscript cross_validation_cv1.R

# Rscript cross_validation_cv2.R

# Rscript cross_validation_cv00.R

## Parent-offspring validation
# Rscript parent_offspring_validation.R

## Sample script
Rscript prediction_validation_sample.R

