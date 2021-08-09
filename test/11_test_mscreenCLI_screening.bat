#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
@REM module load Anaconda3/5.3.0



@REM echo "##############################################"
@REM echo "TESTING SCREENING+PREPARE VINA"
@REM python mscreen.py  screen -l ../data/ligands -r ../data/receptor -c ../data/config_vina_ex2.txt -o ../data/out-vina-screening-prepare -d vina -v -p 
@REM echo SCREENING+PREPARE VINA OK

@REM echo "##############################################"
@REM echo TESTING SCREENING+PREPARE PLANTS
@REM python mscreen.py screen  -l ../data/ligands -r ../data/receptor -c ../data/test_plants_speed1.txt -o ../data/out-plants-screening-prepare -v -d plants -p
@REM echo"SCREENING+PREPARE PLANTS

@REM echo "##############################################"
@REM echo TESTING SCREENING+PREPARE AUTODOCK
@REM python mscreen.py screen  -l ../data/ligands -r ../data/receptor -c ../data/test_autodock.txt -o ../data/out-autodock-screening-prepare -v -d autodock -p
@REM echo PREPARING SCREENING+PREPARE AUTODOCK

@REM echo "##############################################"
@REM echo TESTING SCREENING+PREPARE AUTODOCKZN
@REM python mscreen.py screen  -l ../data/ligands -r ../data/receptor -c ../data/test_autodockzn.txt -o ../data/out-autodockzn-screening-prepare -v -d autodockzn -p
@REM echo PREPARING SCREENING+PREPARE AUTODOCKZN


echo "##############################################"
echo "----------------------------------------------"
echo "----------------------------------------------"
echo "##############################################"

echo "##############################################"
echo "TESTING SCREENING VINA"
python mscreen.py  screen -l ../data/prepare/prepared_ligands_vina -r ../data/prepare/prepared_receptors_vina -c ../data/config_vina_ex2.txt -o ../data/out-vina-screening -d vina -v  
echo SCREENING VINA OK

echo "##############################################"
echo TESTING SCREENING PLANTS
python mscreen.py screen  -l ../data/prepare/prepared_ligands_plants -r ../data/prepare/prepared_receptors_plants -c ../data/config_plants_speed4.txt -o ../data/out-plants-screening -v -d plants 
echo"SCREENING PLANTS

echo "##############################################"
echo TESTING SCREENING AUTODOCK
python mscreen.py screen  -l ../data/prepare/prepared_ligands_autodock -r ../data/prepare/prepared_receptors_autodock -c ../data/config_vina_ex2.txt -o ../data/out-autodock-screening -v -d autodock 
echo PREPARING SCREENING AUTODOCK

echo "##############################################"
echo TESTING SCREENING AUTODOCKZN
python mscreen.py screen  -l ../data/prepare/prepared_ligands_autodock -r ../data/prepare/prepared_receptors_autodock_zn -c ../data/config_vina_ex2.txt -o ../data/out-autodockzn-screening -v -d autodockzn 
echo PREPARING SCREENING AUTODOCKZN

