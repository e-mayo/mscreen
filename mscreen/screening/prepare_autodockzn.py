
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 03:31:00 2021

AutodockZn preparation routines

@author: e-mayo
"""



import os, sys, re
from pathlib import  Path
from shutil import which

sys.path.append('autodocktools_prepare_py3k')
sys.path.append('mscreen/autodocktools_prepare_py3k')

from AutoDockTools.Utilities24.prepare4 import prepare_ligand4, prepare_receptor4, prepare_dpf42, prepare_gpf4


from screening.prepare_autodock import fix_gpf4_file
from screening.prepare_autodock import prepare_autodock42_parameter_file
from screening.prepare_autodock import run_autogrid
from screening.prepare_autodock import run_autodock
from screening.prepare_autodock import prepare_receptor_autodock

from screening.prepare_gpf4zn import prepare_gpf4zn
from screening.zinc_pseudo import zinc_pseudo


def prepare_grid_parameter_file(receptor_filename, ligand_filename,output_gpf_filename,
                   list_parameters_str=None,**kwargs):
    """
    This function call prepare_gpf4zn routine from prepare_gpf4zn.py.
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
    parameter_file_path = Path("screening/AD4Zn.dat")
    if not list_parameters_str:
        list_parameters_str=[f"parameter_file={parameter_file_path.absolute()}"]
        # list_parameters_str=[f'parameter_file="{parameter_file_path}"']
        
    prepare_gpf4zn(receptor_filename = receptor_filename, ligand_filename = ligand_filename, output_gpf_filename = output_gpf_filename, list_parameters_str=list_parameters_str, **kwargs)
    fix_gpf4_file(output_gpf_filename,receptor_filename)


    
def prepare_receptor_autodockzn(receptor, outputfilename=None, **kwargs):
    if not outputfilename:
        outputfilename = Path(receptor.parent / (str(receptor.stem) + '.pdbqt'))
    receptor_prep = prepare_receptor_autodock(receptor,
                                              outputfilename=outputfilename,
                                              **kwargs)
    # receptor_zn = outputfilename.parent / f"{receptor.stem}_zn.pdbqt"
    receptor_zn = zinc_pseudo(receptor_prep, outputfilename)
    return receptor_zn
    