echo '##############################################'
echo TESTING PREPARING VINA
python mscreen.py prepare -bk vina -l ../data/ligands -r ../data/receptor -o ../data/prepare -v
echo PREPARING VINA OK
echo "##############################################"
echo TESTING PREPARING PLANTS
python mscreen.py prepare -bk plants -l ../data/ligands -r ../data/receptor -o ../data/prepare -v
echo PREPARING VINA PLANTS
@REM # python mscreen.py prepare -bk vina -l ../data/PfA-M1_2a/ligands -r ../data/PfA-M1_2a/receptor -o ../data/PfA-M1_2a/prepare 
@REM # python mscreen.py prepare -bk plants -l ../data/PfA-M1_2a/ligands -r ../data/PfA-M1_2a/receptor -o ../data/PfA-M1_2a/prepare
