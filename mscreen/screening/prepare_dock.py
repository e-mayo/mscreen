import sys, os
# sys.path.append('mscreen')


import numpy as np
from pathlib import Path
from shutil import which

sys.path.append(str(Path(sys.path[0]).parent)) #uncomment when working in this file
from screening.pdb_utils import extract_near_residues_as_pdb
import screening.lig_unique_name as lig_unique_name




def dockprep_chimera(mol_file, mol_out=None):
    if not mol_out:
        mol_out = f"{mol_file.split('.')[0]}_prep.mol2"
    DOCKPREP_HOME = Path(__file__).parent / 'dockprep.py'
    command = f"chimera.exe --nogui {mol_file} 'dockprep.py' {mol_out}"
    os.system(command)
    

def prepare_receptor_dock(rec_file, lig_file, sph_selector_cut_off = 10.0, ref=None, out_folder = None):
        """
        Prepare receptor for DOCK6. The receptor must be preprocessed 
        (delete solvent, add hiydrogens, fill missing residues/atoms). 
        Prepare means:
        - generate receptor surface  [dms]
        - calculate spheres  [sphgen] 
        - select sphere within x Angstrom from the binding site [sphere_selector]
        - box generation [showbox]
        - grid generation [grid]
        """
        rec_file = Path(rec_file)
        rec_name = rec_file.stem
        lig_file = Path(lig_file)

        if not out_folder:
            out_folder = rec_file.parent
        else:
            out_folder = Path(out_folder)

        if not ref:
            ref = ref_center = calc_mol_center(lig_file)
        
        # sometimes dock file when atoms id go beyond 9999 or residues id go beyond 999
        # this extract a proteing fragment [frag_name]
        # within radii of the  referece coordinates 
        frag_name = f"frag-{rec_name}"
        frag_file = out_folder / f"{frag_name}.pdb"
        radii = 20

        extract_near_residues_as_pdb(rec_file, frag_file, ref , radii)
        
        print('calculating surface')
        dms_file = out_folder / f"{frag_name}.dms"
        run_dms(frag_file, dms_file)
        
        print('generating spheres')
        sph_file = out_folder / f"{frag_name}.sph"
        # write_INSPH(dms_file, sph_file) 
        # run_spheres(sph_file, dms_file) # use run_sphgen_cpp insted
        run_sphgen_cpp(dms_file, sph_file)
        
        print('selecting spheres')
        sel_file = out_folder / f'{sph_file.stem}-sele.sph'
        run_spheres_selector(sph_file, lig_file, sel_file, cut_off = sph_selector_cut_off)
        

        print("writing box")
        showbox_in = out_folder / f'{rec_file.stem}-showbox.in'
        box_file = out_folder / f'{rec_file.stem}-box'
        write_boxin(sel_file, frag_file,box_file=box_file, showbox_file=showbox_in)
        run_showbox(showbox_in=showbox_in)

        print("writing grid")
        grid_in = out_folder / f'{rec_file.stem}-grid.in'
        grid_out = out_folder / f'{rec_file.stem}-grid.out'
        grid_prefix = out_folder / f'{rec_file.stem}-grid'
        write_gridin(rec_file, box_file, grid_file_name=grid_in, grid_prefix=grid_prefix)
        run_grid(grid_in, grid_out)
        # box_file.unlink()
        return True

def run_dms(rec_file, dms_file):
    """
    Run dms surface generation
    rec_file:: receptor file name (rec.pdb)
    dms_file:: receptor file name (rec.dms)
    """
    rec_file = Path(rec_file)
    rec_name = rec_file.stem
    os.system(f"dms {rec_file}  -a -g {dms_file}.log -n -w 1.4 -v -o {dms_file}")

def run_spheres(sph_file, INSPH=None):
    cwd = Path.cwd()
    if INSPH:
        os.chdir(INSPH.parent)
    try:
        os.remove(sph_file)
    except FileNotFoundError:
        pass
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
    if INSPH:
        os.chdir(cwd)

def run_sphgen_cpp(dms_file, sph_file):
    # sphgen_cpp = './docking_executable/linux/sphegen_cpp/sphgen_cpp'
    clean_dms_file(dms_file)
    command = f'sphgen_cpp -i {dms_file} -o {sph_file}'   
    os.system(command)
    

def get_dock_env_var():
    DOCK_EXE  = Path(which("dock6"))
    DOCK_HOME = DOCK_EXE.parents[1]
    PARAMETERS_HOME = DOCK_HOME / 'parameters'
    return DOCK_EXE, DOCK_HOME, PARAMETERS_HOME 
    
def write_gridin(rec_file, box_file, grid_file_name="grid.in", grid_prefix='grid'):
    DOCK_EXE, DOCK_HOME, PARAMETERS_HOME = get_dock_env_var()
    with open(f"{grid_file_name}", "w") as f:
        lines = ["compute_grids                             yes\n",
                "grid_spacing                              0.4\n",
                "output_molecule                           no\n",
                "contact_score                             no\n",
                "energy_score                              yes\n",
                "energy_cutoff_distance                    999\n",
                "atom_model                                a\n",
                "attractive_exponent                       6\n",
                "repulsive_exponent                        9\n",
                "distance_dielectric                       yes\n",
                "dielectric_factor                         4\n",
                "bump_filter                               yes\n",
                "bump_overlap                              0.75\n",
                f"receptor_file                             {rec_file}\n",
                f"box_file                                  {box_file}\n",
                f"vdw_definition_file {PARAMETERS_HOME / 'vdw_AMBER_parm99.defn'}\n",
                # "vdw_definition_file                      {DOCK_HOME}/parameters/vdw_AMBER_parm99.defn\n",
                f"score_grid_prefix                        {grid_prefix}\n"]
        f.writelines(lines)
        
def write_boxin(sph_file, rec_file, box_size=4, box_file=None, showbox_file="showbox.in"):
    if not box_file:
        rec_file = Path(rec_file)
        rec_name = rec_file.stem
        box_file = rec_file.parent / f"{rec_name}.box.pdb"
    with open(f"{showbox_file}", "w") as f:
        lines = ["Y\n",           # -> Says yes, we want the box to be generated                                
                 f"{box_size}\n", # -> How far in angstroms should the box appear from our selected spheres                          
                 f"{sph_file}\n", # -> Location of your our selected spheres file
                 "1\n",
                 f"{box_file}"]
        f.writelines(lines)


def write_INSPH(dms_file, sph_file=None, s_clas=0.0,
                 max_radii=4.0,
                 min_radii=1.4):
    """Write the INSPH file
    dms_file:: receptor dms file (rec.dms)
    sph_file:: sphere file name default = {rec_name}.sph
    s_clas:: Steric clashing default: 0.0
    max_radii:: Max sphere radii default: 4.0
    min_radii:: Min sphere radii default: 0.0
    """
    dms_file = Path(dms_file)
    dms_name = dms_file.stem
    INSPH = dms_file.parent / f"{dms_name}-INSPH"
    
    if not sph_file:
        sph_file = dms_file.parent / f'{dms_file}.sph'

    
    with open(INSPH, "w") as f:
        lines = [f"{dms_file}\n",
                "R\n",
                "X\n",
                f"{s_clas}\n",
                f"{max_radii}\n",
                f"{min_radii}\n",
                f"{sph_file}"]
        f.writelines(lines)


def write_showspheres(sph_file):
    sph_file = Path(sph_file)
    sph_file = sph_file.stem
    rec_name = sph_file.split('.')[0]
    with open("showsphere.in", "w") as f:
        lines = [f"{sph_file}\n",
                "-1\n",
                "N\n",
                "clustertemp\n",
                "N\n"]
        f.writelines(lines)

def run_grid(grid_file, grid_out):
    os.system(f"grid -i {grid_file} -o {grid_out}")

def run_spheres_selector(sph_file, lig_file, sel_file=None, cut_off=10.0):
    if not sel_file:
        sel_file = sph_file.parent / f'{sph_file.stem}-sele.sph'
    os.system(f"sphere_selector {sph_file} {lig_file} {cut_off}")
    os.rename('selected_spheres.sph', sel_file)

def clean_dms_file(sph_file):
    with open(sph_file, 'r') as f:
        text = f.read()
    text = text.replace("*", "")
    with open(sph_file, 'w') as f:
        f.write(text)

        

def run_showbox(showbox_in="showbox.in"):
    os.system(f"showbox < {showbox_in}")
    
def calc_mol_center(mol2_file):
    """
    Calculate molecular center.

    Parameters
    ----------
    mol2_file : str
        Reference molecule.

    Returns
    -------
    center : np.array [x,y,z]
        Center of reference calculated as the average of all atoms coordinates.

    """
    with  open(mol2_file, "r") as f:
        lines = f.readlines()
    atom_xyz_list = []
    
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            if "ATOM" in line:
                atom_record = True
                continue
            else:
                atom_record = False
        if atom_record:
            atom_line = line.split()
            
            atom_id  = int(atom_line[0])
            atom_name = str(atom_line[1])
            coord_x  = float(atom_line[2])
            coord_y  = float(atom_line[3])
            coord_z  = float(atom_line[4])
            elem     = str(atom_line[5])
            chain_id = int(atom_line[6])
            res_name = str(atom_line[7])
            charge   = float(atom_line[8])
            atom_xyz_list.append([coord_x, coord_y, coord_z])
    atoms = np.array(atom_xyz_list)
    return atoms.mean(axis=0)

def write_preprocess_ligand(lig_file,file_name='dock.lig.in',out_prefix='lig'):
    """
    Write a dock in file for preprocessing the ligand.
    ref:https://github.com/rizzolab/DOCK6_Screening_Protocols/blob/master/run.001.lig_clean_am1bcc.csh

    Parameters
    ----------
    lig_file : str
        Ligand mol2 file.
    file_name : str, optional
        Dock.in file name. The default is 'dock.lig.in'.

    Returns
    -------
    None.

    """
    
    DOCK_EXE, DOCK_HOME, PARAMETERS_HOME = get_dock_env_var()
    
    lines = ['conformer_search_type                                       rigid\n',
            'use_internal_energy                                          no\n',
            f'ligand_atom_file                                             {lig_file}\n',
            'limit_max_ligands                                            no\n',
            'skip_molecule                                                no\n',
            'read_mol_solvation                                           no\n',
            'calculate_rmsd                                               no\n',
            'use_database_filter                                          no\n',
            'orient_ligand                                                no\n',
            'bump_filter                                                  no\n',
            'score_molecules                                              no\n',
            f'vdw_defn_file                                               {PARAMETERS_HOME /"vdw_AMBER_parm99.defn"}\n',
            f'flex_defn_file                                              {PARAMETERS_HOME /"flex.defn"}\n',
            f'flex_drive_file                                             {PARAMETERS_HOME /"flex_drive.tbl"}\n',
            f'ligand_outfile_prefix                                       {out_prefix}\n',
            'write_orientations                                           no\n',
            'num_scored_conformers                                        1\n',
            'rank_ligands                                                 no\n']
    
    
    with open(f'{file_name}', 'w') as f:
        for i in lines:
            print(i)
            f.write(i)
        # f.writelines(lines)
        
def run_dock_6(in_file,out_file):
    os.system(f'dock6 -i {in_file} -o {out_file}')        

def dockprep_ligand(lig_file, lig_out=None):
    if not lig_out:
        lig_out = f"{lig_file.split('.')[0]}_prep.mol2"
    command = f"chimera --nogui {lig_file} dockprep.py {lig_out}"
    os.system(command)

def prepare_ligand_dock(lig_file,out_file=None):
    
    lig_file = Path(lig_file)
    lig_name = lig_file.stem
    
    
    if not out_file:
        prepared_file = lig_file.parent / f'{lig_name}-prep.mol2'
        preprocess_ligand = lig_file.parent / f'{lig_name}-unique_resn.mol2'
    else:
        out_file = Path(out_file)
        preprocess_ligand = out_file.parent / f'{lig_name}-unique_resn.mol2'
        prepared_file = out_file

    lig_unique_name.unique_residue(lig_file, preprocess_ligand)        
    dockprep_chimera(preprocess_ligand, prepared_file)
    preprocess_ligand.unlink()
    

    
if __name__ == "__main__":
    from time import time
    
    
    lig_file = "../data/ligands/lig_3ebh.mol2"
    out_file = "../data/prepare/rec_3ebh.mol2"
    # prepare_ligand_dock(lig_file, out_file)
    
    rec_file = "../data/receptor/rec_3ebh.pdb"
    rec_name = rec_file.split()[0]
    prep_rec_name = f'{rec_name}-prep.pdb'
    t0 = time()
    dockprep_chimera(rec_file, prep_rec_name)
    print(f'dockprep rec preparation time {time()-t0}')
    
    t0 = time()
    out_folder="../data/prepare/"
    prepare_receptor_dock(rec_file, lig_file, sph_selector_cut_off=10.0, out_folder=out_folder)
    print(f'prepare_receptor_dock rec preparation time {time()-t0}')
    
    # ref_center = calc_mol_center(lig_file)

    # print("writing grid")
    # rec_name = rec_file.split('.')[0]
    # frag_name = f"frag-{rec_name}"
    # frag_file = f"{frag_name}.pdb"
    # box_file = f"{frag_name}.box.pdb"
    # write_gridin(rec_file, box_file, grid_file_name='grid.in',grid_prefix="new_grid__.grid")
    # grid_file = "grid.in"
    # grid_out = "grid.out"
    # run_grid(grid_file, grid_out)