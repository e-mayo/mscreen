# -*- coding: utf-8 -*-
"""
Created on Wed May  5 18:17:43 2021

@author: o_o
"""

import os
import re
from pdb_utils import extract_near_residues_as_pdb, renumbering
# prepare_receptor require that AmberTools is installed
# usually DOCK6 came with AmberTools



def run_dms(rec_file, dms_file):
    """
    Run dms surface generation
    rec_file:: receptor file name (rec.pdb)
    dms_file:: receptor file name (rec.dms)
    """
    os.system(f"dms rec_file  -a -g {rec_name}.dms.log -n -w 1.4 -v -o dms_file")


def write_INSPH(dms_file, sph_file=None, s_clas=0.0,
                 m_radii=4.0,
                 min_radii=1.4):
    """Write the INSPH file
    dms_file:: receptor dms file (rec.dms)
    sph_file:: sphere file name default = {rec_name}.sph
    s_clas:: Steric clashing default: 0.0
    max_radii:: Max sphere radii default: 4.0
    min_radii:: Min sphere radii default: 0.0
    """
    if not sph_file:
        sph_file = {rec_name}.sph
    with open("INSPH", "w") as f:
        lines = [f"{dms_file}\n",
                "R\n",
                "X\n",
                f"{s_clas}\n",
                f"{max_radii}\n",
                f"{min_radii}\n",
                f"{sph_file}"]
        f.writelines(lines)

def write_showspheres(sph_file):
    with open("showsphere.in", "w") as f:
        lines = [f"{rec_name}.sph\n",
                "-1\n",
                "N\n",
                "clustertemp\n",
                "N\n"]
        f.writelines(lines)

def write_boxin(sel_sph_file, rec_file, box_size=4, box_file=None, showbox_file="showbox.in"):
    if not box_file:
        rec_name = rec_file.split('.')[0]
        box_file = f"{rec_name}.box.pdb"
    with open(f"{showbox_file}", "w") as f:
        lines = ["Y\n",# -> Says yes, we want the box to be generated                                
                 f"{box_size}\n", # -> How far in angstroms should the box appear from our selected spheres                          
                 f"{sph_file}\n", # -> Location of your our selected spheres file
                 "1\n",
                 f"{box_file}"]
        f.writelines(lines)

def write_gridin(rec_file, box_file, grid_file_name="grid.in"):
    with open(f"{grid_file_name}", "w") as f:
        lines = ["compute_grids                             yes\n",
                "grid_spacing                              0.4\n",
                "output_molecule                           no\n",
                "contact_score                             no\n",
                "energy_score                              yes\n",
                "energy_cutoff_distance                    9999\n",
                "atom_model                                a\n",
                "attractive_exponent                       6\n",
                "repulsive_exponent                        9\n",
                "distance_dielectric                       yes\n",
                "dielectric_factor                         4\n",
                "bump_filter                               yes\n",
                "bump_overlap                              0.75\n",
                f"receptor_file                             {rec_file}\n",
                f"box_file                                  {box_file}\n",
                f"vdw_definition_file /mnt/f/eduardo/programas/docking/DOCK6/dock6/parameters/vdw_AMBER_parm99.defn"
                # "vdw_definition_file                      /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/vdw_AMBER_parm99.defn\n",
                "score_grid_prefix                        grid\n"]

def run_grid(grid_file, grid_out):
    os.system(f"grid -i {grid_file} -o {grid_out}")

def run_spheres_selector(sph_file, lig_file, cut_off=10.0):
    os.system(f"sphere_selector {sph_file} {lig_file} {cut_off}")

def run_showbox(showbox_in="showbox.in"):
    os.system(f"showbox < {showbox_in}")

def run_spheres():
    try:
        os.remove('OUTSPH')
    except FileNotFoundError:
        pass
    try:
        os.remove('temp1.ms')
    except FileNotFoundError:
        pass
    try:
        os.remove('temp3.atc')
    except FileNotFoundError:
        pass 
    os.system("sphgen")


def run_box():
    os.system("showbox < box.in")


# def run_grid():
#     os.system("grid -i grid.in -o grid.out")

def dms2xyz(dms_file,xyz_file):
    pattern = re.compile('\w{3}\s*\d+\w\s*\w+\s*([-\d\.]*)\s*([-\d\.]*)\s*([-\d\.]*).*')
    with open(dms_file, 'r') as f:
        text = f.read()
    textn = pattern.sub(r"H \1 \2 \3", text)
    with open(xyz_file, 'w') as f:
        f.write(textn)

in_pdb_file = 'frag_3ebh_noH.pdb'
# in_pdb_file = 'frag_3ebh_noH.pdb'

def renumbering_residues(in_pdb_file,out_pdb_file='python.pdb'):
    with open(in_pdb_file, 'r') as f:
        text = f.read()
    
    text = re.sub('ANISOU.*\n', '', text)
    pattern = re.compile('(ATOM\s+\d+\s+\w+\s+\w{3}\s\w\s*)(\d*)(.*)')
    # pattern = re.compile('ATOM\s+\d+\s+\w+\s+\w{3}\s\w\s*(\d*)')
    m = pattern.findall(text)
    first_resid = m[0][1]
    a_id = 2
    for a in m:
        if first_resid != a[1]:
            a_id+=1
            first_resid = a[1]
            print(f'a_id:{a_id}--first_resid:{first_resid}--a[1]:{a[1]}')
        text = text.replace(f'{a[0]}{a[1]}{a[2]}', f'{a[0]}{a_id}{a[2]}')
        # if a_id == 20:break
    
    # print(text[:2000])
    with open(out_pdb_file, 'w') as f:
        f.write(text)

# def run_dock_6.7(num):
#     os.system("dock6 -i anchor_and_grow_"+str(num)+".in -o anchor_and_grow_"+str(num)+".out")

# renumbering_residues('frag_3ebh_noH.pdb','frag_3ebh_noHr.pdb')
rec_file = "frag_3ebh_noHr.pdb"
dms_file = "frag_3ebh_noHr.dms"
sph_file = "frag_3ebh_noHr.sph"
lig_file = "lig_3ebh.mol2"
box_file = "frag_3ebh_noHr.sph"

extract_near_residues_as_pdb(rec_file,f'sele-{rec_file}',[10.15, 10.08, 4.02 ],20)
chains = process_pdb_file(pdb_file)
protein = chains[0]
# chains = process_pdb_file(f'{pdb_file.split(".")[0]}-rn.pdb')
# protein = chains[0]
selection = select([10.15, 10.08, 4.02 ], 20, protein, byres=True)
# atoms =  
write_pdb('selection_byres.pdb', selection)
chains = process_pdb_file('selection_byres.pdb')
protein = chains[0]
protein = renumbering(protein,resi=True,atomi=True,save_as=f'selection_byres-rn.pdb')
write_pdb('selection_byres_rn.pdb', selection)

# run_dms(rec_file, dms_file)
# write_INSPH(dms_file,sph_file)
# run_spheres()

# run_spheres_selector(sph_file, lig_file, cut_off=10.0)
sel_file = "selected_spheres.sph"

# write_boxin(sel_file, rec_file)

# run_showbox()

write_gridin(rec_file, box_file)
grid_file = "grid.in"
grid_out = "grid.out"
# run_grid(grid_file, grid_out)


# dms2xyz('frag_3ebh_noH.dms','frag_3ebh_noH.xyz')
# import os
# os.removedirs('003.spheres',)
# os.mkdir('003.spheres')
# os.chdir('003.spheres')


