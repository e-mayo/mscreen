# mscreen
Mscreen try to be a common user interface for different ligand-protein docking programs.

# Usage
It's quite simple to use: prepare a ligand and receptor for screening

## Ligand and receptor preparation
python mscreen.py prepare -bk [backend] -l [ligand_folder] -r [receptor_folder] -o [otput_folder]

## Running virtual screening
python mscreen.py screen -bk plants -l [ligand_folder] -r [receptor_folder] -o [otput_folder] -c [conf_file] -log [log_file]

python mscreen.py analysis pending

# Ependencies
Make sure you have the following dependencies:
rdkit 
numpy



