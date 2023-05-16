#!/bin/sh

#SBATCH --job-name=JOBNAME
#SBATCH --output=LOGFILENAME
#SBATCH --error=ERRORFILENAME
#SBATCH --partition=htc
#SBATCH --ntasks=NTASKS
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000
#SBATCH --time=02:00:00

JOBS
