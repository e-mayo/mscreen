#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
@REM module load Anaconda3/5.3.0

@REM echo '##############################################'
@REM echo TESTING ANALYSIS VINA
python mscreen.py analysis -bk vina   -vs ../data/out-test-vina_ex1  -t short
python mscreen.py analysis -bk vina   -vs ../data/out-test-vina_ex1  -t full
@REM echo ANALYSIS VINA OK

@REM echo PREPARING VINA OK
@REM echo "##############################################"
@REM echo TESTING ANALYSIS PLANTS
python mscreen.py analysis -bk plants -vs ../data/out-test-plants_speed4 -t short
python mscreen.py analysis -bk plants -vs ../data/out-test-plants_speed4 -t full
@REM echo ANALYSIS PLANTS OK



