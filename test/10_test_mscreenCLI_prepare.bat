echo "##############################################"
echo "TESTING PREPARING VINA"
python mscreen.py  prepare -l ../data/ligands -r ../data/receptors -o ../data/prepare -d vina -v
echo PREPARING VINA OK
echo "##############################################"
echo TESTING PREPARING PLANTS
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d plants
echo "PREPARING VINA PLANTS"
echo "##############################################"
echo TESTING PREPARING AUTODOCK
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodock
echo "PREPARING AUTODOCK"
echo "##############################################"
echo TESTING PREPARING AUTODOCK ZN
python mscreen.py prepare  -l ../data/ligands -r ../data/receptors -o ../data/prepare -v -d autodockzn
echo "PREPARING AUTODOCK ZN"

