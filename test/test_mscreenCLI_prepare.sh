#$ -S /bin/sh -cwd
# $  -pe mpi.singlehost 12
export OMP_NUM_THREADS=14
module load Anaconda3/5.3.0

echo '##############################################'
echo TESTING PREPARING VINA
python mscreen.py prepare -bk vina -l ../data/ligands -r ../data/receptor -o ../data/prepare -v
echo PREPARING VINA OK
echo "##############################################"
echo TESTING PREPARING PLANTS
python mscreen.py prepare -bk plants -l ../data/ligands -r ../data/receptor -o ../data/prepare -v
echo PREPARING VINA PLANTS
# python mscreen.py prepare -bk vina -l ../data/PfA-M1_2a/ligands -r ../data/PfA-M1_2a/receptor -o ../data/PfA-M1_2a/prepare 
# python mscreen.py prepare -bk plants -l ../data/PfA-M1_2a/ligands -r ../data/PfA-M1_2a/receptor -o ../data/PfA-M1_2a/prepare
