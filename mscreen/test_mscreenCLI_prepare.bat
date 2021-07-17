echo "##############################################"
@REM echo "TESTING PREPARING VINA"
@REM python mscreen.py  prepare -l ../data/ligands -r ../data/receptor -o ../data/prepare -d vina -v
@REM echo PREPARING VINA OK
echo "##############################################"
echo TESTING PREPARING AUTODOCK
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodock
echo "PREPARING AUTODOCK"