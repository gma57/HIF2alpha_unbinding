#!/bin/bash
#SBATCH --job-name=Holo_p2s4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=72:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --mail-type=END
#SBATCH --mail-user=ris130@pitt.edu

module purge
module load gcc/8.2.0 openmpi/4.0.3
#module load amber/20_cuda-11.1.1
module load amber/22-beta

#pmemd.cuda -O -i prod.in -o prod.out -p BOUND.prmtop -c 4_eq3.rst -r prod.rst -x prod.nc -inf prod.nfo

#pmemd.cuda -O -i prod.in -o prod2.out -p BOUND.prmtop -c prod.rst -r prod2.rst -x prod2.nc -inf prod2.nfo

#pmemd.cuda -O -i prod.in -o prod3.out -p BOUND.prmtop -c prod2.rst -r prod3.rst -x prod3.nc -inf prod3.nfo

pmemd.cuda -O -i prod.in -o prod4.out -p BOUND.prmtop -c prod3.rst -r prod4.rst -x prod4.nc -inf prod4.nfo

