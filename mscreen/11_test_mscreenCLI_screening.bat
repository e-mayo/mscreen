#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
@REM module load Anaconda3/5.3.0



@REM @echo ##############################################
@REM @echo "TESTING SCREENING+PREPARE VINA"
@REM python mscreen.py  screen -l ../data/ligands -r ../data/receptors -c ../data/config_vina_ex2.txt -o ../data/out-vina-screening-prepare -d vina -v -p 
@REM @echo SCREENING+PREPARE VINA OK

@echo ##############################################
@echo TESTING SCREENING+PREPARE PLANTS
python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_plants_speed1.txt -o ../data/out-plants-screening-prepare -v -d plants -p
@echo "SCREENING+PREPARE PLANTS

@echo ##############################################
@echo TESTING SCREENING+PREPARE AUTODOCK
python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_vina_ex1.txt -o ../data/out-autodock-screening-prepare -v -d autodock -p
@echo PREPARING SCREENING+PREPARE AUTODOCK

@echo ##############################################
@echo TESTING SCREENING+PREPARE AUTODOCKZN
python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_vina_ex1.txt -o ../data/out-autodockzn-screening-prepare -v -d autodockzn -p
@echo PREPARING SCREENING+PREPARE AUTODOCKZN


@echo ##############################################
@echo ----------------------------------------------
@echo ----------------------------------------------
@echo ##############################################

@echo ##############################################
@echo "TESTING SCREENING VINA"
python mscreen.py  screen -l ../data/prepare/prepare_ligands_vina -r ../data/prepare/prepare_receptors_vina -c ../data/config_vina_ex2.txt -o ../data/out-vina-screening -d vina -v  
@echo SCREENING VINA OK

@echo ##############################################
@echo TESTING SCREENING PLANTS
python mscreen.py screen  -l ../data/prepare/prepare_ligands_plants -r ../data/prepare/prepare_receptors_plants -c ../data/test_plants_speed1.txt -o ../data/out-vina-plants-screening -v -d plants 
@echo SCREENING PLANTS

@echo ##############################################
@echo TESTING SCREENING AUTODOCK
python mscreen.py screen  -l ../data/prepare/prepare_ligands_autodock -r ../data/prepare/prepare_receptors_autodock -c ../data/test_autodock.txt -o ../data/out-vina-autodock-screening -v -d autodock 
@echo PREPARING SCREENING AUTODOCK

@echo ##############################################
@echo TESTING SCREENING AUTODOCKZN
python mscreen.py screen  -l ../data/prepare/prepare_ligands_autodock -r ../data/prepare/prepare_receptors_autodock -c ../data/test_autodockzn.txt -o ../data/out-vina-autodockzn-screening -v -d autodockzn 
@echo PREPARING SCREENING AUTODOCKZN

