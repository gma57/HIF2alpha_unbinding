#!/bin/bash
#SBATCH --job-name=eq4apo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=72:00:00
#SBATCH --cluster=invest
#SBATCH --partition=lchong
##SBATCH --mail-type=END,BEGIN,FAIL
##SBATCH --mail-user=gma57@pitt.edu


module purge
module load gcc/8.2.0 openmpi/4.0.3 amber/20_cuda-11.1.1

TOP=bound_new.prmtop
CRD=bound_new.inpcrd


echo $AMBERHOME

echo "First step: minimizing the whole system"
mkdir 1_min 2_eq 3_eq 4_eq
cd 1_min
pmemd.cuda -O -i ../min1.in -p ../${TOP} -c ../${CRD} -o 1_min.out -r 1_min.rst -x 1_min.nc -inf 1_min.nfo
cd ../
echo "Done first step"


echo "Step 2: heating up the system"
cd 2_eq
pmemd.cuda -O -i ../2_eq1.in -p ../${TOP} -c ../1_min/1_min.rst -o 2_eq1.out -r 2_eq1.rst -x 2_eq1.nc -inf 2_eq1.nfo -ref ../1_min/1_min.rst
cd ../
echo "Done with heating the system"


echo "Step 3: pressurizing the system"
cd 3_eq
pmemd.cuda -O -i ../3_eq2.in -p ../${TOP} -c ../2_eq/2_eq1.rst -o 3_eq2.out -r 3_eq2.rst -x 3_eq2.nc -inf 3_eq2.nfo -ref ../2_eq/2_eq1.rst 
cd ../
echo "Done step 3"



echo "Step 4: equilibrating whole system"
cd 4_eq
pmemd.cuda -O -i ../4_eq3.in -p ../${TOP} -c ../3_eq/3_eq2.rst -o 4_eq3.out -r 4_eq3.rst -x 4_eq3.nc -inf 4_eq3.nfo
echo "Done step 4"
