# mscreen
Mscreen try to be a common user interface for different ligand-protein docking programs.

# Installation
1- Install python 3.7+
2- Install rdkit
3- Install numpy
git clone https://github.com/e-mayo/mscreen.git


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

## Docking program supported
* ### Autodock Vina 
Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455–461. https://doi.org/10.1002/jcc.21334
* ### FWAVina
Li J, Song Y, Li F, Zhang H, Liu W. FWAVina: A novel optimization algorithm for     protein-ligand docking based on the fireworks algorithm. Comput Biol Chem. 2020 Oct;88:107363. doi: 10.1016/j.compbiolchem.2020.107363. Epub 2020 Aug 20. PMID: 32861160.
* ### GWOVina
Wong, KM, Tai, HK, Siu, SWI. GWOVina: A grey wolf optimization approach to rigid and flexible receptor docking. Chem Biol Drug Des. 2021; 97: 97– 110. https://doi.org/10.1111/cbdd.13764

* ### PSOVina
Ng MC, Fong S, Siu SW. PSOVina: The hybrid particle swarm optimization algorithm for protein-ligand docking. J Bioinform Comput Biol. 2015 Jun;13(3):1541007. doi: 10.1142/S0219720015410073. Epub 2015 Feb 10. PMID: 25800162.
* ## QVina-2.1
“Fast, Accurate, and Reliable Molecular Docking with QuickVina 2” Amr Alhossary, Stephanus Daniel Handoko, Yuguang Mu, and Chee-Keong Kwoh. Bioinformatics (2015) 31 (13): 2214-2216. DOI:10.1093/bioinformatics/btv082
* ### QVina-w
“Protein-Ligand Blind Docking Using QuickVina-W With Inter-Process Spatio-Temporal Integration” Nafisa M. Hassan, Amr A. Alhossary, Yuguang Mu & Chee-Keong Kwoh. Nature Scientific Reports 7(1) (2017). DOI:10.1038/s41598-017-15571-7
* ### PLANTS
Korb, O.; Stützle, T.; Exner, T. E. PLANTS: Application of Ant Colony Optimization to Structure-Based Drug Design, Lecture Notes in Computer Science 4150, 247-258 (2006)
# Who we are?
We are SOAN from the Laboratorio de Sintesis Organica at The University of Havana

