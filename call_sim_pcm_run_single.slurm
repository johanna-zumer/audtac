#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --qos bbdefault
#SBATCH --job-name sim_pcm
#SBATCH --mem 2G # 4GB
#SBATcH --mem-per-cpu 2G
#SBATCH --time 48:00:00
#SBATCH --mail-type=NONE  #ALL
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --account=zumerj01

# Because --export=NONE is set, we need to unset
# SLURM_EXPORT_ENV, otherwise multi-node MPI jobs will fail

unset SLURM_EXPORT_ENV
set -e
module purge; module load bluebear
module load apps/matlab/r2017a

matlab -nodisplay -r "sim2_pcm_run($tr,$sm)"

