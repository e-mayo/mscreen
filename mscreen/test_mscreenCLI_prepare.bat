echo "##############################################"
echo "TESTING PREPARING VINA"
python mscreen.py  prepare -l ../data/ligands -r ../data/receptor -o ../data/prepare -d vina -v
echo PREPARING VINA OK
echo "##############################################"
echo TESTING PREPARING PLANTS
python mscreen.py prepare  -l ../data/ligands -r ../data/receptor -o ../data/prepare -v -d plants
echo "PREPARING VINA PLANTS"