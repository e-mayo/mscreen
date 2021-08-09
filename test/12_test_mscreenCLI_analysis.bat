#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
@REM module load Anaconda3/5.3.0

@REM echo '##############################################'
echo TESTING ANALYSIS VINA
python mscreen.py analysis -d vina   -i ../data/out-vina-screening  -t short
python mscreen.py analysis -d vina   -i ../data/out-vina-screening  -t full
echo ANALYSIS VINA OK

echo "##############################################"
echo TESTING ANALYSIS AUTODOCK
python mscreen.py analysis -d autodock -i ../data/out-autodock-screening -t short
python mscreen.py analysis -d autodock -i ../data/out-autodock-screening -t full
echo ANALYSIS AUTODOCK OK

echo "##############################################"
echo TESTING ANALYSIS AUTODOCKZN
python mscreen.py analysis -d autodockzn -i ../data/out-autodockzn-screening -t short
python mscreen.py analysis -d autodockzn -i ../data/out-autodockzn-screening -t full
echo ANALYSIS AUTODOCKZN OK


echo "##############################################"
echo TESTING ANALYSIS PLANTS
python mscreen.py analysis -d plants -i ../data/out-plants-screening -t short
python mscreen.py analysis -d plants -i ../data/out-plants-screening -t full
echo ANALYSIS PLANTS OK

echo "##############################################"
echo TESTING ANALYSIS DOCK
python mscreen.py analysis -d dock -i ../data/out-dock-screening -t short
python mscreen.py analysis -d dock -i ../data/out-dock-screening -t full
echo ANALYSIS DOCK OK
echo "##############################################"


