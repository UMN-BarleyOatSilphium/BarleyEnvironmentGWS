#!/bin/bash

#PBS -l walltime=12:00:00,mem=24gb,nodes=1:ppn=1
#PBS -N historical_ec_timeframe_selection
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/Scripts/EnvironmentalCovariables/

module load R/3.5.2_mkl

# Run the historical covariate timeframe selection
Rscript covariate_phenotypic_modeling_location.R

