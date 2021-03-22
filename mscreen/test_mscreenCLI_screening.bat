#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
@REM module load Anaconda3/5.3.0

@REM echo '##############################################'
@REM echo TESTING SCREEN VINA
@REM python mscreen.py screen -l ../data/prepare/prepared_ligands_vina -r ../data/prepare/prepared_receptors_vina -o ../data/out-test-vina_ex1  -c ../data/config_vina_ex8.txt  -log ../data/workshop_vina_ex8.log -v -d vina

@REM echo PREPARING VINA OK
echo "##############################################"
echo TESTING SCREEN PLANTS
python mscreen.py screen -l ../data/prepare/prepared_ligands_plants -r ../data/prepare/prepared_receptors_plants -o ../data/out-test-plants_speed4 -c ../data/config_plants_speed4.txt -log ../data/test_plants_speed4.log -d plants -v 
echo PREPARING PLANTS OK



