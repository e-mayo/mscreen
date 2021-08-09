#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8
module load Anaconda3/5.3.0

echo '##############################################'
echo TESTING SCREEN VINA
python mscreen.py screen -d vina   -l ../data/prepare/prepared_ligands_vina \
                                    -r ../data/prepare/prepared_receptors_vina \
                                    -o ../data/out-test-vina_ex1 \
                                    -c ../data/config_vina_ex8.txt \
                                    -log ../data/PfA-M1_2a_vina_ex8.log\
                                    -v

echo PREPARING VINA OK
echo "##############################################"
echo TESTING SCREEN dock
python mscreen.py screen -d dock -l ../data/prepare/prepared_ligands_dock -r ../data/prepare/prepared_receptors_dock -o ../data/out-test-dock_speed4  -c ../data/config_dock_flex_sample.txt -log ../data/test_dock_flex.log -v
echo PREPARING dock OK


