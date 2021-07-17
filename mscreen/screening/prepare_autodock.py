
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 17:35:34 2020

Prepare ligand from python 3 using the python 2.7 interpreter of ADFRsuite

@author: emayo
"""
import os, sys, re
from pathlib import  Path
from shutil import which

sys.path.append('autodocktools_prepare_py3k')
sys.path.append('mscreen/autodocktools_prepare_py3k')
from screening.prepare_vina import *
from AutoDockTools.Utilities24.prepare4 import prepare_ligand4, prepare_receptor4, prepare_dpf42, prepare_gpf4



def prepare_autodock42_parameter_file(receptor_filename, ligand_filename,dpf_filename=None, **kwargs):
    """
    This function call the python wrapper in AutodockTools/Utils24/prepare4 for prepare_dpf42.py script.
    It prepare the dpf file for autodock4.
    """
    receptor_filename = str(receptor_filename)
    ligand_filename = str(ligand_filename)
    if dpf_filename:
        dpf_filename = str(dpf_filename)
    prepare_dpf42(receptor_filename=receptor_filename, ligand_filename=ligand_filename,dpf_filename=dpf_filename, **kwargs)
    fix_dpf_file(dpf_filename, ligand_filename)

def fix_dpf_file(dpf_filename,lig):
    dpf_filename= Path(dpf_filename)
    lig = Path(lig)
    with open(dpf_filename,'r') as f:
        dpf_lines = f.readlines()
    pattern = re.compile('^(\w+)\s(.+?)\s+(.*)')
    
    parent = str(dpf_filename.parent).replace('\\','/')
    lig = str(lig).replace('\\','/')
    new_dpf_lines = []
    for line in dpf_lines:
        if line.startswith("move "):
            line = pattern.sub(rf"\1 {lig} \3", line)
        if line.startswith("fld "):
            line= pattern.sub(rf"\1 {parent}/\2 \3", line)
        if line.startswith("map "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)
        if line.startswith("elecmap "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)
        if line.startswith("desolvmap "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)

        new_dpf_lines.append(line)
        
    with open(dpf_filename,'w+') as f:
        f.writelines(new_dpf_lines)
    return dpf_filename    

def prepare_grid_parameter_file(receptor_filename, ligand_filename,output_gpf_filename, **kwargs):
    """
    This function call the python wrapper in AutodockTools/Utils24/prepare4 for prepare_gpf4.py script.
    It prepare the gpf file for autogrid4.
    Parameters
    ----------
    receptor_filename : str, optional
        DESCRIPTION. The default is None.
    ligand_filename : str, optional
        ligand_filename (.pdbq format). The default is None.
    list_filename : TYPE, optional
        DESCRIPTION. The default is None.
    gpf_filename : str, optional
        reference_gpf_filename. The default is None.
    output_gpf_filename : str, optional
        output_gpf_filename. The default is None.
    flexres_filename : TYPE, optional
        flexres_filename. The default is None.
    directory : TYPE, optional
        directory of ligands to use to set types. The default is None.
    list_parameters_str : list, optional
        parameter=newvalue. For example: ["ligand_types='HD,Br,A,C,OA'","ligand_types='HD,Br,A,C,OA'"] or ["gridcenter='2.5,6.5,-7.5'""]. The default is [].
    verbose : TYPE, optional
        DESCRIPTION. The default is None.
    center_on_ligand : TYPE, optional
        boolean to center grids on center of ligand. The default is False.
    size_box_to_include_ligand : TYPE, optional
        boolean to NOT size_box_to_include_ligand. The default is True.
    npts_increment : int, optional
        increment npts in all 3 dimensions by this integer. The default is 0.
    ligand_types_defined : TYPE, optional
        DESCRIPTION. The default is False.
    verbose : TYPE, optional
        print stuff. The default is False.

    Returns
    -------
    None.

    """
    receptor_filename = str(receptor_filename)
    ligand_filename = str(ligand_filename)
    if output_gpf_filename:
        output_gpf_filename = str(output_gpf_filename)
    prepare_gpf4(receptor_filename=receptor_filename,
                ligand_filename=ligand_filename,
                output_gpf_filename=output_gpf_filename,**kwargs) 
    fix_gpf4_file(output_gpf_filename,receptor_filename)

def fix_gpf4_file(gpf_filename,rec):
    gpf_filename= Path(gpf_filename)
    rec = Path(rec)
    with open(gpf_filename,'r') as f:
        gpf_lines = f.readlines()
    pattern = re.compile('^(\w+)\s(.+?)\s+(.*)')
    
    parent = str(gpf_filename.parent).replace('\\','/')
    rec = str(rec).replace('\\','/')
    new_gpf_lines = []
    for line in gpf_lines:
        if line.startswith("receptor "):
            line= pattern.sub(rf"\1 {rec} \3", line)
        if line.startswith("gridfld "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)
        if line.startswith("map "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)
        if line.startswith("elecmap "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)
        if line.startswith("dsolvmap "):
            line = pattern.sub(rf"\1 {parent}/\2 \3", line)

        new_gpf_lines.append(line)
        
    with open(gpf_filename,'w+') as f:
        f.writelines(new_gpf_lines)
    return gpf_filename
            
   

        
    
    
    
def prepare_receptor_autodock(receptor, outputfilename=None, **kwargs):
    if not outputfilename:
        outputfilename = Path(receptor.parent / (str(receptor.stem) + '.pdbqt'))
    receptor = str(receptor)
    try:
        prepare_receptor4(receptor,outputfilename=outputfilename, **kwargs)
    except:
        print(f'fail at receptor {outputfilename.name}')
        return None
    return outputfilename

def prepare_ligand_autodock(ligand, outputfilename=None, **kwargs):
    if not outputfilename:
        outputfilename = Path(ligand.parent / (str(ligand.stem) + '.pdbqt'))
    ligand = str(ligand)
    
    try:
        prepare_ligand4(ligand,outputfilename=outputfilename, **kwargs)
    except:
        print(f'fail at ligand {outputfilename.name}')
        return None
    return outputfilename
    
def get_autodock_env_var():
    AUTOGRID_EXE = Path(which("autogrid4"))
    AUTODOCK_EXE = Path(which("autodock4"))
    AUTODOCK_HOME = AUTODOCK_EXE.parents[1]
    return AUTOGRID_EXE, AUTODOCK_EXE

def run_autogrid(parameter_filename, log_filename=None):
    
    command = f"autogrid4 -p {parameter_filename}"
    """
    # command = f'"{AUTOGRID_EXE}" -p {parameter_filename}'
    windows path has spaces this cause a problem when the executable
    is placed on a folder with spaces in its path. one way to solve this is using ""
    """
    
    if log_filename:
        command += f" -l {log_filename}"
    os.system(command)

def run_autodock(parameter_filename, log_filename=None,
                keep_original_residue=False,
                ignore_header_checking=False):
    parameter_filename = str(parameter_filename)
    
    command = f'autodock4 -p "{parameter_filename}" '
    if log_filename:
        log_filename = str(log_filename)
        command += f"-l {log_filename} "
    if keep_original_residue:
        command += "-k "
    if ignore_header_checking:
        command += "-i "
    os.system(command)

if __name__ == "__main__":
    #AutdoDock42 docking test:
    from pathlib import Path
    ligand = Path("data/ligands/lig_5y19.mol2")
    out_lig = Path("data/prepare/prepared_ligands_autodock") / (str(ligand.stem) + '.pdbqt')
    receptor = Path("data/receptors/rec_5y19.pdb")
    out_rec  = Path("data/prepare/prepared_receptors_autodock") / (str(receptor.stem) + '.pdbqt')
    
    # prepare_ligand_autodock(ligand, outputfilename=out_lig)
    # prepare_receptor_autodock(receptor, outputfilename=out_rec)

    outputPath = Path("data/prepare") 
    output_gpf_filename = outputPath / "test.gpf"
    prepare_grid_parameter_file(receptor_filename=out_rec,ligand_filename=out_lig, output_gpf_filename=output_gpf_filename)
    run_autogrid(parameter_filename=output_gpf_filename, log_filename=None)
    
    dpf_filename = outputPath / f"{out_rec.stem}-{out_lig.stem}.dpf"
    prepare_autodock42_parameter_file(receptor_filename=out_rec, # abug mix lig_name and rec_name
                                      ligand_filename=out_lig,
                                      dpf_filename=dpf_filename)
    # automatically create  lig_5y19_rec_5y19.dlg logfile
    run_autodock(dpf_filename, log_filename=None,
                            keep_original_residue=False,
                            ignore_header_checking=False)
    