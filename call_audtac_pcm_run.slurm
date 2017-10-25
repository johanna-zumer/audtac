#!/bin/bash

# Example multi-job script designed to submit a single job script to 
# the scheduler multiple times, each time using the index of a different 
# EEG subject.
#
# This attemps to use the cluster to run repetitive jobs
# (such as running the same analysis on different subject) in parallel. 
#
# In this example, the script that runs the analysis for each subject 
# uses 1 variable:
# 	SUBID = the subject id number

# Loop through subjects and submit each one to the scheduler:
trind=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25) 
# smind=(1 2 3 4) 
for t in ${trind[@]}; do
	#for s in ${smind[@]}; do
#		sbatch --export=tr=${t},sm=${s} call_sim_pcm_run_single.slurm
		sbatch --export=tr=${t} call_audtac_pcm_run_single.slurm
	#done
done

