#!/bin/sh

#SBATCH --job-name=JOBNAME
#SBATCH --output=LOGFILENAME
#SBATCH --error=ERRORFILENAME
#SBATCH --partition=htc
#SBATCH --ntasks=NTASKS
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000
#SBATCH --time=23:00:00
#SBATCH --licenses=sps

start=$(date +%s)

JOBS

end=$(date +%s)
echo "Time JOBNAME: $((($end-$start)/60)) mins"
