#!/bin/bash
#SBATCH --job-name=lpath
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=11:59:59
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --mail-user=mls333@pitt.edu
#SBATCH --mail-type=BEGIN,END

source env.sh || exit 1

## Step 1: extract trajectories
lpath extract -we -W ../multi.h5 -A ../ANALYSIS/OVERALL/assign.h5 -ss 0 -ts 1 -p --stride 1 --auxdata='dist_lig' -rw --trace-basis

## Step 2: match trajectories and clusterize them
lpath match -we -ra reassign_custom.reassign_custom --remove-ends --condense --output-pickle reassigned.pickle -dP -1 -pdt 30 -dh --debug -ex -W ../multi.h5

#lpath match -we -ra reassign_custom.reassign_custom --remove-ends --condense --output-pickle reassigned.pickle -dN -dF succ_traj/distmat.npy
