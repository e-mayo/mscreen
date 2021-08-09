echo "##############################################"
@REM echo "TESTING PREPARING VINA"
@REM python mscreen.py  prepare -l ../data/ligands -r ../data/receptor -o ../data/prepare -d vina -v
@REM echo PREPARING VINA OK
@REM echo "##############################################"
@REM echo TESTING PREPARING AUTODOCK
@REM python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodock
@REM echo "PREPARING AUTODOCK"
echo "##############################################"
echo TESTING PREPARING AUTODOCK ZN
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodockzn
echo "PREPARING AUTODOCK ZN"