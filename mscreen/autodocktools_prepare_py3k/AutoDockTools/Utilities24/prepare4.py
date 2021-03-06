import os, sys 

from MolKit import Read

import MolKit.molecule
import MolKit.protein
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
from AutoDockTools.MoleculePreparation import AD4LigandPreparation


def prepare_ligand4(ligand_filename,
                    verbose = None,
                    add_bonds = False,
                    repairs = "",
                    charges_to_add = None,
                    preserve_charge_types='',
                    cleanup  = "nphs_lps",
                    allowed_bonds = "backbone",
                    root = 'auto',
                    outputfilename = None,
                    check_for_fragments = False,
                    bonds_to_inactivate = "",
                    inactivate_all_torsions = False,
                    attach_nonbonded_fragments = False,
                    attach_singletons = False,
                    assign_unique_names = False,
                    mode = 'automatic',
                    dictionary= None):
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
    
    if ligand_filename:
        ligand_filename = ligand_filename
        if verbose: print('set ligand_filename to ', ligand_filename)
        ligand_ext = os.path.splitext(os.path.basename(ligand_filename))[1] # eg .mol2  
    if verbose:
        if verbose: print('set verbose to ', True)
    if outputfilename:
        if verbose: print('set outputfilename to ', outputfilename)
    if dictionary:
        if verbose: print('set dictionaryto ', dictionary)
    if repairs:
        if verbose: print('set repairs to ', repairs)
    if not charges_to_add:
        charges_to_add = 'gasteiger'
        if verbose: print('Adding gasteiger')
    else: 
        if verbose: print('Do not add charges gasteiger')
    if preserve_charge_types:
        if ligand_ext in ['.mol2']:
            if preserve_charge_types =='':
                if preserve_charge_types.find(',') == -1:
                    preserve_charge_types = preserve_charge_types
                else:
                    preserve_charge_types = []
                    for subT in preserve_charge_types.split(','):
                        preserve_charge_types.append(subT)
        if verbose: print('preserve initial charges on ', preserve_charge_types)
    if cleanup:
        if verbose: print('set cleanup to merge ', cleanup)
    if allowed_bonds:
        if verbose: print('allow ', allowed_bonds, 'bonds set to rotate')
    if root:
        if verbose: print('set root to ', root)
    if check_for_fragments:
        if verbose: print('set check_for_fragments to True')
    if mode:
        if verbose: print('set mode to ', mode)
    if bonds_to_inactivate:
        if verbose: print('set bonds_to_inactivate to ', bonds_to_inactivate)
    if inactivate_all_torsions:
        if verbose: print('set inactivate_all_torsions to ', inactivate_all_torsions)
    if attach_nonbonded_fragments:
        if verbose: print('set attach_nonbonded_fragments to ', attach_nonbonded_fragments)
    if attach_singletons:
        if verbose: print('set attach_singletons to ', attach_singletons)
    if assign_unique_names:
        if verbose: print('set assign_unique_names to ', assign_unique_names, ' newname is original name plus its index(1-based')
    
    if not ligand_filename:
        raise ValueError('prepare_ligand4: ligand filename must be specified.')
        return None
    if not outputfilename:
        out_name = os.path.splitext(os.path.basename(ligand_filename))[0] + '.pdbqt'
        outputfilename = os.path.join(os.path.dirname(ligand_filename), os.path.basename(out_name))
    if attach_singletons:
        attach_nonbonded_fragments = True
        if verbose: print("using attach_singletons so attach_nonbonded_fragments also")
    
    mols = Read(ligand_filename)
    if verbose: print('read ', ligand_filename)
    mol = mols[0]
    if len(mols)>1:
        if verbose: 
            print("more than one molecule in file")
        #use the one molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose:
                    print("mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms")
    coord_dict = {}
    for a in mol.allAtoms: coord_dict[a] = a.coords
    if assign_unique_names:  # added to simplify setting up covalent dockings 8/2014
        for at in mol.allAtoms:
            if mol.allAtoms.get(at.name) >1:
                at.name = at.name + str(at._uniqIndex +1)
        if verbose:
            print("renamed %d atoms: each newname is the original name of the atom plus its (1-based) uniqIndex" %(len(mol.allAtoms)))

    mol.buildBondsByDistance()
    if charges_to_add is not None:
        preserved = {}
        preserved_types = preserve_charge_types.split(',') 
        for t in preserved_types:
            if not len(t): continue
            try:
                ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
                for a in ats:
                    if a.chargeSet is not None:
                        preserved[a] = [a.chargeSet, a.charge]
            except AttributeError:
                ats = mol.allAtoms.get(lambda x: x.element==t)
                for a in ats:
                    if a.chargeSet is not None:
                        preserved[a] = [a.chargeSet, a.charge]
            print(" preserved = ", end=' ') 
            for key, val in list(preserved.items()):
                print("key=", key)
                print("val =", val)
           

    if verbose:
        print("setting up LPO with mode=", mode, end=' ')
        print("and outputfilename= ", outputfilename)
        print("and check_for_fragments=", check_for_fragments)
        print("and bonds_to_inactivate=", bonds_to_inactivate)
    LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add, 
                        cleanup, allowed_bonds, root, 
                        outputfilename=outputfilename,
                        dict=dictionary, check_for_fragments=check_for_fragments,
                        bonds_to_inactivate=bonds_to_inactivate, 
                        inactivate_all_torsions=inactivate_all_torsions,
                        attach_nonbonded_fragments=attach_nonbonded_fragments,
                        attach_singletons=attach_singletons)
    #do something about atoms with too many bonds (?)
    #FIX THIS: could be peptide ligand (???)
    #          ??use isPeptide to decide chargeSet??
    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in list(preserved.items()):
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]
            if verbose: print("set charge on ", atom.full_name(), " to ", atom.charge)
    if verbose: print("returning ", mol.returnCode) 
    bad_list = []
    for a in mol.allAtoms:
        if a in list(coord_dict.keys()) and a.coords!=coord_dict[a]: 
            bad_list.append(a)
    if len(bad_list):
        print(len(bad_list), ' atom coordinates changed!')    
        for a in bad_list:
            print(a.name, ":", coord_dict[a], ' -> ', a.coords)
    else:
        if verbose: print("No change in atomic coordinates")
    if mol.returnCode!=0: 
        sys.stderr.write(mol.returnMsg+"\n")
    return(mol.returnCode)
    
    
    
    

def prepare_receptor4(receptor_filename,
                        repairs = '',
                        charges_to_add = 'gasteiger',
                        preserve_charge_types=None,
                        cleanup  = "nphs_lps_waters_nonstdres",
                        outputfilename = None,
                        mode = 'automatic',
                        delete_single_nonstd_residues = None,
                        dictionary = None,
                        unique_atom_names = False,
                        verbose = None):
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
    ### Initialize the function parameters as in the prepare_receptor4.py script
    if receptor_filename:
        if verbose: print('set receptor_filename to ', receptor_filename)
    if verbose:
        verbose = True
        if verbose: print('set verbose to ', True)
    if outputfilename:
        if verbose: print('set outputfilename to ', outputfilename)
    if repairs:
        if verbose: print('set repairs to ', repairs)
    if not charges_to_add:
        if verbose: print(f'Adding {charges_to_add}')
    # if not charges_to_add:
        # charges_to_add = 'gasteiger'
        # if verbose: print('Adding gasteiger')
    if preserve_charge_types:   
        preserve_charge_types = str(preserve_charge_types)
        import re    
        pattern = r'[\'\"\[\]]'
        r = re.compile(pattern)
        preserve_charge_types = r.sub('', preserve_charge_types)
        if verbose: print('preserve initial charges on ', preserve_charge_types)
    if cleanup:
        if verbose: print('set cleanup to ', cleanup)
    if delete_single_nonstd_residues:
        delete_single_nonstd_residues = True
        if verbose: print('set delete_single_nonstd_residues to True')
    if mode:
        if verbose: print('set mode to ', mode)
    if dictionary:
        if verbose: print('set dictionary to ', dictionary)
    if unique_atom_names:
        unique_atom_names = True
        if verbose: print('set unique_atom_names to ', unique_atom_names)

    if not receptor_filename:
        raise ValueError('prepare_receptor4: receptor filename must be specified.')
        return None
    
    if not outputfilename:
        out_name = os.path.splitext(os.path.basename(receptor_filename))[0] + '.pdbqt'
        outputfilename = os.path.join(os.path.dirname(receptor_filename), os.path.basename(out_name))
    #what about nucleic acids???

    mols = Read(receptor_filename)
    if verbose: print('read ', receptor_filename)
    mol = mols[0]
    if unique_atom_names:  # added to simplify setting up covalent dockings 8/2014
        for at in mol.allAtoms:
            if mol.allAtoms.get(at.name) >1:
                at.name = at.name + str(at._uniqIndex +1)
        if verbose:
            print("renamed %d atoms: each newname is the original name of the atom plus its (1-based) uniqIndex" %(len(mol.allAtoms)))        
    preserved = {}
    has_autodock_element = False
    print('#'*72)
    print(f'charges_to_add is not None {charges_to_add is not None}')
    print(f'charges_to_add: {charges_to_add}')
    if charges_to_add is not None and preserve_charge_types is not None:      
        if hasattr(mol, 'allAtoms') and not hasattr(mol.allAtoms[0], 'autodock_element'):
            file_name, file_ext = os.path.splitext(receptor_filename)
            if file_ext == '.pdbqt':
                has_autodock_element = True
        if preserve_charge_types is not None and has_autodock_element==False:
            print('prepare_receptor4: input format does not have autodock_element SO unable to preserve charges on ' + preserve_charge_types)
            print('exiting...')
            return None  
        preserved_types = preserve_charge_types.split(',') 
        if verbose: print("preserved_types=", preserved_types)
        for t in preserved_types:
            if verbose: print('preserving charges on type->', t)
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            if verbose: print("preserving charges on ", ats.name)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]

    if len(mols)>1:
        if verbose: print("more than one molecule in file")
        #use the molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose: print("mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms")
    mol.buildBondsByDistance()
    alt_loc_ats = mol.allAtoms.get(lambda x: "@" in x.name)
    len_alt_loc_ats = len(alt_loc_ats)
    if len_alt_loc_ats:
        print("WARNING!", mol.name, "has",len_alt_loc_ats, ' alternate location atoms!\nUse prepare_pdb_split_alt_confs.py to create pdb files containing a single conformation.\n')

    if verbose:
        print("setting up RPO with mode=", mode, end=' ')
        print("and outputfilename= ", outputfilename)
        print("charges_to_add=", charges_to_add)
        print("delete_single_nonstd_residues=", delete_single_nonstd_residues)

    RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, 
                        cleanup, outputfilename=outputfilename,
                        preserved=preserved, 
                        delete_single_nonstd_residues=delete_single_nonstd_residues,
                        dict=dictionary)    

    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in list(preserved.items()):
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]
            

if __name__ == '__main__':
    from pathlib import Path
    receptor = 'C:\\Users\\o_o\\Documents\\research\\postgrad\\SOAN-UH\\E.Coli_P.Falsiparum\\structures\\03_docking_tutorial\\data\\receptor_3b2x.pdb'
    ligand = 'C:\\Users\\o_o\\Documents\\research\\postgrad\\SOAN-UH\\E.Coli_P.Falsiparum\\structures\\03_docking_tutorial\\data\\lig_3b2x.mol2'
    # outputfilename = 'C:\\Users\\o_o\\Documents\\research\\postgrad\\SOAN-UH\\E.Coli_P.Falsiparum\\structures\\03_docking_tutorial\\data\\lig_3b2x.pdbqt'
    prepare_receptor4(receptor, verbose=True)
    # prepare_ligand4(ligand,verbose=True,charges_to_add=None)#, outputfilename=outputfilename)
    