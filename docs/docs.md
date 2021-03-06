#Mscree documentation
Hello finally mscreen is out:
mscreen is an open source program for running virtual screening. It help you running virtual screening using different open source programs. We choose the huge family of autodock vina. 
Docking programs avalaibles:
    vina 
    gwovina:
    Wong, K. M., Tai, H. K., & Siu, S. W. (2021). GWOVina: A grey wolf optimization approach to rigid and flexible receptor docking. Chemical Biology & Drug Design, 97(1), 97-110.
    psovina2.0: 
    Hio Kuan Tai, Siti Azma Jusoh, and Shirley W. I. Siu*. Chaos-embedded particle swarm optimization approach for protein-ligand docking and virtual screening. Journal of Cheminformatics  10:62, 2018.    
    fwavina:
    Jin Li, Yongping Song, Fajin Li, Henggui Zhang, Weichao Liu, FWAVina: A novel optimization algorithm for protein-ligand docking based on the fireworks algorithm, Computational Biology and Chemistry,  88,2020
    qvina-2.1:
    Alhossary, A., Handoko, S. D., Mu, Y., & Kwoh, C. K. (2015). Fast, accurate, and reliable molecular docking with QuickVina 2. Bioinformatics, 31(13), 2214-2216.
    qvina-w:
    Hassan, N. M., Alhossary, A. A., Mu, Y., & Kwoh, C. K. (2017). Protein-ligand blind docking using QuickVina-W with inter-process spatio-temporal integration. Scientific reports, 7(1), 1-13.
    SMINA:
    Koes, D. R., Baumgartner, M. P., & Camacho, C. J. (2013). Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise. Journal of chemical information and modeling, 53(8), 1893-1904.

#What's our goal?
Make open source virtual screening as easy as possible.

How does mscreen work?
If you want run a virtual screening run the follow comand.
python mscreen.py screen -l path_to_ligand_folder
                        -r path_to_receptor_folder
                        -o path_to_the_output_folder
                        -c path_to_the_configuration_fila

#Tutorial
Setting mscreen up.
Download mscreen link
Install dependecies.
Add mscreen.py to your path.

#Running a virtual screening.

##setting up the screening workplace
Create the following directory tree.

virtual_screening-|- ligands
                  |- receptors
                  |- conf.txt

Minimal configuration file must contain the following parameters.
center_x =  14.2
center_y = 43.6
center_z = 36.9
size_x = 31.3
size_y = 20.4
size_z = 21.5
exhaustiveness = 8
num_modes = 10

##screening a library of ligand against one or more receptor
Run the following command and voila!
mscreen screen -l ligand_folder [The path to ligands Default: ligands]
                -r receptor_folder [The path to receptors Default: receptors]
                -o output_folder [The path to the output folder Default: out]
                -c configuration_fila [The name of the configuration file for docking']
                -bk vina [The docking software to use, choose one of the next 'fwavina', 'gwovina',     'ledock', 'psovina2', 'qvina2.1', 'qvina-w', 'smina', 'vina' default: qvina-w]

After running it you will find into output folder the following files and folder.
One folder for each receptor used and the configuration file used.
Inside the receptor folder you will find one folder for each ligand docked.
If you want use a deferent docking program add the -bk flag and the name of the program. Available programs 'fwavina', 'gwovina', 'ledock', 'psovina', 'qvina2', 'qvina-w', 'smina', 'vina'. Unfortunatelly for windows only vina is available. 

##analysis of the virtual screening result.
mscreen analysis -a full/short [The type of analysis default: short]
                 -o output_folder [The name of the output folder. Default: out]
                 -rn [Report file name default:results.txt]

###short analysis
output a single file inside each receptor folder with the ligands names and score sorted from low to high

###full analysis
run a short analysis but also output a single file sdf file containing all structures docked with the following information:
ligand name, ligand score, ligand cluster, is representative


