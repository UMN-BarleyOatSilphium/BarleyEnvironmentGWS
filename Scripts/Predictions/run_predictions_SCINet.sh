#!/bin/bash

#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -t 12:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem-per-cpu 8G   # maximum memory per node
#SBATCH -p short # GPU node with 21 day maximum, but can be killed anytime
#SBATCH --job-name="barley_environment_GWS_predictions"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Change the working directory
cd /project/gifvl_vaccinium/barley_work/BarleyEnvironmentGWS/

module load r-intel/3.6.1

# Leave-one-environment-out
Rscript Scripts/Predictions/01_environmentl_predictions_varcomp_SCINet.R
