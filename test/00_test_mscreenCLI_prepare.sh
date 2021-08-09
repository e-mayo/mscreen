#$ -S /bin/sh -cwd
# $  -pe mpi.singlehost 12
# export OMP_NUM_THREADS=14
# module load Anaconda3/5.3.0

echo "##############################################"
echo "TESTING PREPARING VINA"
python mscreen.py  prepare -l ../data/ligands -r ../data/receptors -o ../data/prepare -d vina -v
echo PREPARING VINA OK

echo "##############################################"
echo TESTING PREPARING PLANTS
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d plants
echo "PREPARING VINA PLANTS"

echo "##############################################"
echo TESTING PREPARING DOCK
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d dock
echo "PREPARING VINA DOCK"

echo "##############################################"
echo TESTING PREPARING AUTODOCK
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodock
echo PREPARING PREPARING AUTODOCK

echo "##############################################"
echo TESTING PREPARING AUTODOCKZN
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodockzn
echo PREPARING PREPARING AUTODOCKZN

