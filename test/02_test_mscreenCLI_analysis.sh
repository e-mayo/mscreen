#$ -S /bin/sh -cwd
# $  -pe mpi.singlehost 12
export OMP_NUM_THREADS=14
module load Anaconda3/5.3.0

echo '##############################################'
echo TESTING ANALYSIS VINA
python mscreen.py analysis -d vina   -i ../data/out-vina-screening-prepare  -t short
python mscreen.py analysis -d vina   -i ../data/out-vina-screening-prepare  -t full
echo ANALYSIS VINA OK


echo "##############################################"
echo TESTING ANALYSIS PLANTS
python mscreen.py analysis -d plants -i ../data/out-plants-screening-prepare -t short
python mscreen.py analysis -d plants -i ../data/out-plants-screening-prepare -t full
echo ANALYSIS PLANTS OK

echo "##############################################"
echo TESTING ANALYSIS DOCK
    python mscreen.py analysis -d dock -i ../data/out-dock-screening -t short
python mscreen.py analysis -d dock -i ../data/out-dock-screening-prepare -t full
echo ANALYSIS DOCK OK
echo "##############################################"
