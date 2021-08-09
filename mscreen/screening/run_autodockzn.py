import os, sys, re
from pathlib import  Path
from shutil import which


sys.path.append('autodocktools_prepare_py3k')
sys.path.append('mscreen/autodocktools_prepare_py3k')
from screening.prepare_vina import *
from AutoDockTools.Utilities24.prepare4 import prepare_ligand4, prepare_receptor4, prepare_dpf42, prepare_gpf4
from screening.prepare_autodock import *
from screening.zinc_pseudo import zinc_pseudo
from prepare_gpf4zn import prepare_gpf4zn


rec = Path("rec_5y1k.pdb")
lig = Path("lig_5y1k.mol2")

lig_prep = lig.parent / f"{lig.stem}.pdbqt"
rec_prep = rec.parent / f"{rec.stem}.pdbqt"
rec_zn = rec.parent / f"{rec.stem}_zn.pdbqt"
# prepare_receptor4(rec)
# prepare_ligand4(lig)
zinc_pseudo(rec_prep,rec_zn)


output_gpf_filename = rec.parent / f"{rec.stem}-{lig.stem}.gpf"
        
prepare_gpf4zn(ligand_filename=lig_prep,
               receptor_filename=rec_zn,
               output_gpf_filename=output_gpf_filename,
               list_parameters_str=["parameter_file=AD4Zn.dat"])

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

run_autogrid(output_gpf_filename, log_filename=None)
    
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

prepare_autodock42_parameter_file(rec_zn, lig_prep)   
run_autodock("lig_5y1k_rec_5y1k_zn.dpf")
# python prepare_gpf4zn.py -l lig_5y1k.pdbqt -r rec_5y1k_TZ.pdbqt -o rec_5y1k_TZ.gpf  –p parameter_file=AD4Zn.dat
# autogrid4 -p rec_5y1k_TZ.gpf
