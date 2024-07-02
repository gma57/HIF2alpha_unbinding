#runs 5 production runs in separate folders

for x in 1 2 3 4 5; do

mkdir $x
cd $x
sed 's/NUM/$x/g' ../run.slurm > run.slurm
sbatch run.slurm
cd ../
done
