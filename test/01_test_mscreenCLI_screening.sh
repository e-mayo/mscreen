#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
module load Anaconda3/5.3.0


# echo "##############################################"
# echo "TESTING SCREENING+PREPARE VINA"
# python mscreen.py  screen -l ../data/ligands -r ../data/receptors -c ../data/config_vina_ex2.txt -o ../data/out-vina-screening-prepare -d vina -v -p 
# echo SCREENING+PREPARE VINA OK

# echo "##############################################"
# echo TESTING SCREENING+PREPARE PLANTS
# python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_plants_speed1.txt -o ../data/out-plants-screening-prepare -v -d plants -p
# echo"SCREENING+PREPARE PLANTS

# echo "##############################################"
# echo TESTING SCREENING+PREPARE  DOCK
# python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_dock_sample.txt -o ../data/out-dock-screening-prepare -v -d dock -p
# echo SCREENING+PREPARE  DOCK

# echo "##############################################"
# echo TESTING SCREENING+PREPARE AUTODOCK
# python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_vina_ex1.txt -o ../data/out-autodock-screening-prepare -v -d autodock -p
# echo PREPARING SCREENING+PREPARE AUTODOCK

# echo "##############################################"
# echo TESTING SCREENING+PREPARE AUTODOCKZN
# python mscreen.py screen  -l ../data/ligands -r ../data/receptors -c ../data/config_vina_ex1.txt -o ../data/out-autodockzn-screening-prepare -v -d autodockzn -p
# echo PREPARING SCREENING+PREPARE AUTODOCKZN


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
python mscreen.py screen  -l ../data/prepare/prepared_ligands_plants -r ../data/prepare/prepared_receptors_plants -c ../data/config_plants_speed1.txt -o ../data/out-plants-screening -v -d plants 
echo"SCREENING PLANTS

echo "##############################################"
echo TESTING SCREENING DOCK
python mscreen.py screen  -l ../data/prepare/prepared_ligands_dock -r ../data/prepare/prepared_receptors_dock -c ../data/config_dock_sample.txt -o ../data/out-dock-screening -v -d dock 
echo"SCREENING DOCK

echo "##############################################"
echo TESTING SCREENING AUTODOCK
python mscreen.py screen  -l ../data/prepare/prepared_ligands_autodock -r ../data/prepare/prepared_receptors_autodock -c ../data/config_vina_ex1.txt -o ../data/out-autodock-screening -v -d autodock 
echo PREPARING SCREENING AUTODOCK

echo "##############################################"
echo TESTING SCREENING AUTODOCKZN
python mscreen.py screen  -l ../data/prepare/prepared_ligands_autodock -r ../data/prepare/prepared_receptors_autodock -c ../data/config_vina_ex1.txt -o ../data/out-autodockzn-screening -v -d autodockzn 
echo PREPARING SCREENING AUTODOCKZN

