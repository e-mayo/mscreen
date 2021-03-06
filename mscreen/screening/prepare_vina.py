# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 17:35:34 2020

Prepare ligand from python 3 using the python 2.7 interpreter of ADFRsuite

@author: emayo
"""
import os, sys
from pathlib import  Path

# cwd = Path.cwd()
# print(cwd)
# os.chdir(Path('autodocktools_prepare_py3k/'))
sys.path.append('autodocktools_prepare_py3k')
# print(f'####################### { Path.cwd()}')
from AutoDockTools.Utilities24.prepare4 import prepare_ligand4, prepare_receptor4
# os.chdir(cwd)

def prepare_ligand_vina(ligand,outputfilename=None, **kwargs):
    """
    ligand_filename : str, Required
        ligand_filename (.pdb or .mol2 or .pdbq format)
    verbose : str, optional
        Verbose output. The default is None.
    add_bonds : str, optional
        DESCRIPTION. The default is False.
    repairs : str, optional
        Type of repairs to make:\n\t\t bonds_hydrogens, bonds, hydrogens. The default is is to do no repairs.
    charges_to_add : str, optional
        Do not add charges default. The default is to add 'gasteiger' charges.
    preserve_charge_strs : str, optional
        Preserve input charges on an atom type, eg -p Zn. The default is not to preserve charges on any specific atom type.
    cleanup : str, optional
        Cleanup type: nphs_lps, nphs, lps, ''. The default is "nphs_lps".
    allowed_bonds : str, optional
        Type(s) of bonds to allow to rotate. Default sets 'backbone' rotatable and 'amide' + 'guanidinium' non-rotatable.
    root : str, optional
        Index for root. The default is 'auto'.
    outputfilename : str, optional
        Output pdbqt_filename. The default is output filename is ligand_filename_stem + .pdbqt.
    check_for_fragments : str, optional
        Check for and use largest non-bonded fragment. The default is False (not to do this).
    bonds_to_inactivate : str, optional
        String of bonds to inactivate composed of 
        of zero-based atom indices eg 5_13_2_10
            will inactivate atoms[5]-atoms[13] bond
            and atoms[2]-atoms[10] bond. The default is "" (not to inactivate any specific bonds).
    inactivate_all_torsions : str, optional
        Inactivate all active torsions. The default is False (leave all rotatable active except amide and guanidinium).
    attach_nonbonded_fragments : str, optional
        Attach all nonbonded fragments. The default is False.
    attach_singletons : str, optional
        Attach all nonbonded singletons.
        NB: sets attach all nonbonded fragments too. The default is False.
    assign_unique_names : str, optional
        Assign each ligand atom a unique name: newname is original name plus its index(1-based). The default is False.
    mode : str, optional
        mode = 'automatic'|'interactive'. The default is 'automatic'.
    dictionary: str, optional
        Dictionary to write strs list and number of active torsions. The default is None.
    """
    if not outputfilename:
        outputfilename = Path(ligand.parent / (str(ligand.stem) + '.pdbqt'))
    ligand = str(ligand)
    
    try:
        prepare_ligand4(ligand,outputfilename=outputfilename, **kwargs)
    except:
        print(f'fail at receptor {outputfilename.name}')
        return None
    return outputfilename
    
def prepare_receptor_vina(receptor,outputfilename=None, **kwargs):
    """
    

    Parameters
    ----------
    receptor_filename : str
        DESCRIPTION.
    repairs : str, optional
        Repairs to make: add bonds and/or hydrogens or checkhydrogens.
        'bonds_hydrogens': build bonds and add hydrogens
        'bonds': build a single bond from each atom with no bonds to its closest neighbor
        'hydrogens': add hydrogens
        'checkhydrogens': add hydrogens only if there are none already
        'None': do not make any repairs
        The default is 'None'.
    charges_to_add : str, optional
        Preserve all input charges ie do not add new charges
        The default is is addition of 'gasteiger' charges.
    preserve_charge_strs : list, optional
        Preserve charges on specific atom types. The default is None.
        eg: prepare_receptor4('receptor.pdb',preserve_charge_strs=['Zn','Mn'])
    cleanup:  str, optional
        Cleanup by merging nphs_lps, nphs, lps, waters, nonstdres.
        'nphs': merge charges and remove non-polar hydrogens
        'lps': merge charges and remove lone pairs
        'waters': remove water residues
        'nonstdres': remove chains composed entirely of residues of
                types other than the standard 20 amino acids
        'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX
        The default is "nphs_lps_waters_nonstdres"
    outputfilename : str, optional
        Output file name. The default is molecule_name.pdbqt.
    mode : str, optional
        'automatic' or 'interactive'
        Default is 'automatic': outputfile is written with no further user input
        Mode. The default is 'automatic'.
    delete_single_nonstd_residues : str, optional
        Delete every nonstd residue from each chain.
        'True': any residue whose name is not in this list:
                ['CYS','ILE','SER','VAL','GLN','LYS','ASN',
                'PRO','THR','PHE','ALA','HIS','GLY','ASP',
                'LEU', 'ARG', 'TRP', 'GLU', 'TYR','MET',
                'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']
        will be deleted from any chain.
        NB: there are no  nucleic acid residue names at all in the list and no metals.
        The default is False which means not to do this
    dictionary : str, optional
        File to contain receptor summary information
        Dictionary. The default is None.
    unique_atom_names : str, optional
        Assign each receptor atom a unique name: newname is original name plus its index(1-based). The default is False.
    verbose : str, optional
        Verbose. The default is None.

    Example:
    -------    
    >>>prepare_receptor4(pdb_file, outputfilename='otu_file.pdbqt', repairs='checkhydrogens')

    """
    if not outputfilename:
        outputfilename = Path(receptor.parent / (str(receptor.stem) + '.pdbqt'))
    receptor = str(receptor)
    try:
        prepare_receptor4(receptor,outputfilename=outputfilename, **kwargs)
    except:
        print(f'fail at receptor {outputfilename.name}')
        return None
    return outputfilename


if __name__ == '__main__':
    print(__name__)