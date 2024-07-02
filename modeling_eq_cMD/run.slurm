#!/bin/bash
#SBATCH --job-name=tmNUM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=72:00:00
#SBATCH --cluster=invest
#SBATCH --partition=lchong
#SBATCH --mail-type=END
#SBATCH --mail-user=gma57@pitt.edu

module purge
module load gcc/8.2.0 openmpi/4.0.3
module load amber/20_cuda-11.1.1
#module load amber/22-beta

pmemd.cuda -O -i ../prod.in -o prod.out -p ../bound_new.prmtop -c ../4_eq/4_eq3.rst -r prod.rst -x prod.nc -inf prod.nfo

