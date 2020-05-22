#!/bin/bash

#PBS -l walltime=08:00:00,mem=48gb,nodes=1:ppn=8
#PBS -N factorial_regression
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/Scripts/EnvironmentalCovariables/

module load R/3.5.2_mkl

# Factorial regression samples
Rscript covariate_phenotypic_modeling_sample.R


