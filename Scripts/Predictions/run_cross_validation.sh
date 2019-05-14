#!/bin/bash

# #PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=24
#PBS -l walltime=24:00:00,mem=64gb,nodes=1:ppn=8
# #PBS -N cross-validation-cv1
# #PBS -N cross-validation-cv2
# #PBS -N parent_offspring_validation
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

## Parent-offspring validation
Rscript parent_offspring_validation.R
