#!/bin/bash

#PBS -l walltime=36:00:00,mem=62gb,nodes=1:ppn=8
#PBS -N loeo_loyo_predictions
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/Scripts/Predictions/

module load R/3.5.0
# module load R/3.5.2_mkl

# Predictions by environmental rank
Rscript leave_one_out_prediction.R

