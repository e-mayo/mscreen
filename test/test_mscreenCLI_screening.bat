#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
@REM module load Anaconda3/5.3.0

@REM echo '##############################################'
@REM echo TESTING SCREEN VINA
@REM python mscreen.py screen -l ../data/prepare/prepared_ligands_vina -r ../data/prepare/prepared_receptors_vina -o ../data/out-test-vina_ex1  -c ../data/config_vina_ex8.txt  -log ../data/workshop_vina_ex8.log -v -d vina

@REM echo PREPARING VINA OK
echo "##############################################"
@REM echo TESTING SCREEN AUTODOCK
@REM python mscreen.py screen -l ../data/prepare/prepared_ligands_autodock -r ../data/prepare/prepared_receptors_autodock -o ../data/out-test-autodock  -log ../data/test_autodock.log -d autodock -v 
@REM echo SCREEN AUTODOCK OK

@REM echo TESTING FULL SCREEN AUTODOCK
@REM python mscreen.py screen -l ../data/ligands -r ../data/receptors -o ../data/out-test-autodock  -log ../data/test_autodock.log -d autodock -v -p
@REM echo SCREEN AUTODOCK OK

@REM echo TESTING FULL SCREEN AUTODOCK ZN
@REM python mscreen.py screen -l ../data/ligands -r ../data/receptors -o ../data/out-test-autodockzn  -log ../data/test_autodockzn.log -d autodockzn -v -p
@REM echo SCREEN AUTODOCK ZN

echo TESTING FULL SCREEN AUTODOCK ZN
python mscreen.py screen -l ../data/prepare/prepared_ligands_autodock -r ../data/prepare/prepared_receptors_autodock -o ../data/out-test-autodockzn  -log ../data/test_autodockzn.log -d autodockzn -v 
echo SCREEN AUTODOCK ZN