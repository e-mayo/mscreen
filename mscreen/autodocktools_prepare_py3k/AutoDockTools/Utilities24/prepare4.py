import os, sys
from pathlib import  Path
sys.path.append('autodocktools_prepare_py3k')
sys.path.append('../../../../mscreen/autodocktools_prepare_py3k')
sys.path.append('mscreen/autodocktools_prepare_py3k')
from MolKit import Read

import MolKit.molecule
import MolKit.protein
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
from AutoDockTools.MoleculePreparation import AD4LigandPreparation

# add # Use prepare_pdb_split_alt_confs.py to create pdb files containing a single conformation. fix
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
    ligand_filename = Path(ligand_filename)
    if not ligand_filename.exists():
        print(f"{ligand_filename} doesn't exist")
        return False

    
    ligand_filename = str(ligand_filename)
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

import _py2k_string as string
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list4
from AutoDockTools.GridParameters import GridParameter4FileMaker
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper           

def prepare_gpf4(receptor_filename=None,
                ligand_filename = None,
                list_filename = None,
                gpf_filename = None,
                output_gpf_filename = None,
                flexres_filename = None,
                directory = None,
                list_parameters_str = None,
                verbose = None,
                center_on_ligand = False,
                size_box_to_include_ligand = True,
                npts_increment = 0,
                ligand_types_defined=False):   
    """

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
    parameters : str, optional
        parameter=newvalue. For example: -p ligand_types='HD,Br,A,C,OA' or p npts='60,60,66' or gridcenter='2.5,6.5,-7.5']". The default is [].
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
    parameters = []
    
    if ligand_filename:
        if verbose: print('ligand_filename=', ligand_filename)
    if receptor_filename:
        if verbose: print('receptor_filename=', receptor_filename)
    if gpf_filename:
        if verbose: print('reference_gpf_filename=', gpf_filename)
    if flexres_filename:
        if verbose: print('flexres_filename=', flexres_filename)
    if output_gpf_filename:
        if verbose: print('output_gpf_filename=', output_gpf_filename)
    if list_parameters_str:
        for parameters_str in list_parameters_str:
            parameters.append(parameters_str)
            if parameters_str.split('=')[0]=="ligand_types": ligand_types_defined = True
        if verbose: print('parameters=', parameters)
    if directory:
        if verbose: print('directory=', directory)
    if center_on_ligand:
        if verbose: print('set center_on_ligand to ', center_on_ligand)
    if not size_box_to_include_ligand:
        if verbose: print('set size_box_to_include_ligand to ', size_box_to_include_ligand)
    if npts_increment:
        npts_increment = int(npts_increment)
        if verbose: print('set npts_increment to ', npts_increment)
        
    if (not receptor_filename) or (ligand_filename is None and directory is None and ligand_types_defined is False):
        print("prepare_gpf4.py: ligand and receptor filenames")
        print("                    must be specified.")
        return None
    
    gpfm = GridParameter4FileMaker(size_box_to_include_ligand=size_box_to_include_ligand,verbose=verbose)
    if gpf_filename is not None:
        gpfm.read_reference(gpf_filename)
    if ligand_filename is not None:
        gpfm.set_ligand(ligand_filename)
    gpfm.set_receptor(receptor_filename)
    if directory is not None:
        gpfm.set_types_from_directory(directory)
    if flexres_filename is not None:
        flexmol = Read(flexres_filename)[0]
        flexres_types = flexmol.allAtoms.autodock_element
        lig_types = gpfm.gpo['ligand_types']['value'].split()
        all_types = lig_types
        for t in flexres_types:
            if t not in all_types:
                all_types.append(t)
        all_types_string = all_types[0]
        if len(all_types)>1:
            for t in all_types[1:]:
                all_types_string = all_types_string + " " + t
        gpfm.gpo['ligand_types']['value'] = all_types_string
    for param_str in parameters:
        if param_str.find("parameter_file")>-1:
            parameters.append("custom_parameter_file=1")
            break
    for p in parameters:
        key,newvalue = string.split(p, '=')
        if key=='gridcenter' and newvalue.find(',')>-1:
            newvalue = newvalue.split(',')
            newvalue = string.join(newvalue)
        kw = {key:newvalue}
        gpfm.set_grid_parameters(*(), **kw)
    #gpfm.set_grid_parameters(spacing=1.0)
    if center_on_ligand is True:
        gpfm.gpo['gridcenterAuto']['value'] = 0
        cenx,ceny,cenz = gpfm.ligand.getCenter()
        gpfm.gpo['gridcenter']['value'] = "%.3f %.3f %.3f" %(cenx,ceny,cenz)
    if npts_increment:
        orig_npts = gpfm.gpo['npts']['value']  #[40,40,40]
        if verbose: print("before increment npts=", orig_npts)
        for ind in range(3):
            gpfm.gpo['npts']['value'][ind] += npts_increment
        if verbose: print("after increment npts =", gpfm.gpo['npts']['value'])
    gpfm.write_gpf(output_gpf_filename)
    return True

import _py2k_string as string
import os.path
from MolKit import Read
from AutoDockTools.DockingParameters import DockingParameters, genetic_algorithm_list4_2, \
                genetic_algorithm_local_search_list4_2, local_search_list4_2,\
                simulated_annealing_list4_2, epdb_list4_2,\
                DockingParameter42FileMaker

from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
import numpy

def prepare_dpf42(receptor_filename = None,
                    ligand_filename = None,
                    dpf_filename = None,
                    template_filename = None,
                    flexres_filename = None,
                    list_parameters_str = None,
                    pop_seed = False,
                    local_search=False,
                    use_simulated_annealing=False,
                    verbose = None,
                    epdb_output = False,
                    about_root_atoms_only=False):
    
    parameters = []
    parameter_list = genetic_algorithm_local_search_list4_2

    if ligand_filename:
        if verbose: print('ligand_filename =', ligand_filename)
    if receptor_filename:
        if verbose: print('receptor_filename =', receptor_filename)
    if flexres_filename:
        if verbose: print('flexres_filename =', flexres_filename)
    if template_filename:
        if verbose: print('template_filename =', template_filename)
    if dpf_filename:
        if verbose: print('output dpf_filename =', dpf_filename)
    if list_parameters_str:
        for parameters_str in list_parameters_str:
            parameters.append(parameters_str)
        if verbose: print('parameters =', parameters)
    if epdb_output:
        if verbose: print('output epdb file')
        parameter_list = epdb_list4_2
    if parameter_list:
        if verbose: print('parameter_list =', parameter_list)
    if local_search:   #parameter_list_to_write
        local_search = 1
        parameter_list = local_search_list4_2
        if verbose: print('parameter_list =', parameter_list)
    if use_simulated_annealing:   #parameter_list_to_write
        parameter_list = simulated_annealing_list4_2
        if verbose: print('parameter_list =', parameter_list)
    if about_root_atoms_only:   #set about to average of coords of root atoms
        if verbose: print('about_root_atoms_only =', about_root_atoms_only)
    if pop_seed:
        pop_seed = True
        if verbose: print('pop_seed =', pop_seed)
    
    if (not receptor_filename) or (not ligand_filename):
        print("prepare_dpf42.py: ligand and receptor filenames")
        print("                    must be specified.")



    #11/2011: fixing local_search bugs:
    # specifically:
    # 1. quaternion0 0 0 0 0
    # 2. dihe0 0 0 0 0 0 <one per rotatable bond>
    # 3. about == tran0
    # 4. remove tstep  qstep and dstep
    # 5. remove ls_search_freq
    #local_search = local_search_list4_2
    #parameter_list =local_search_list4_2
    dm = DockingParameter42FileMaker(verbose=verbose, pop_seed=pop_seed)
    if template_filename is not None:  #setup values by reading dpf
        dm.dpo.read(template_filename)
    dm.set_ligand(ligand_filename)
    dm.set_receptor(receptor_filename)
    if flexres_filename is not None:
        flexmol = Read(flexres_filename)[0]
        flexres_types = flexmol.allAtoms.autodock_element
        lig_types = dm.dpo['ligand_types']['value'].split()
        all_types = lig_types
        for t in flexres_types:
            if t not in all_types:
                all_types.append(t)
        all_types_string = all_types[0]
        if len(all_types)>1:
            for t in all_types[1:]:
                all_types_string = all_types_string + " " + t
                if verbose: print("adding ", t, " to all_types->", all_types_string)
        dm.dpo['ligand_types']['value'] = all_types_string
        dm.dpo['flexres']['value'] = flexres_filename
        dm.dpo['flexres_flag']['value'] = True
    #dm.set_docking_parameters( ga_num_evals=1750000,ga_pop_size=150, ga_run=20, rmstol=2.0)
    kw = {}
    for p in parameters:
        key,newvalue = string.split(p, '=')
        #detect string reps of lists: eg "[1.,1.,1.]"
        if key=='parameter_file':
            if key in parameter_list:
                print("removing parameter_file keyword")
                parameter_list.remove('parameter_file')
            parameter_list.insert(1, key)
            dm.dpo['custom_parameter_file']['value']=1
        if newvalue[0]=='[':
            nv = []
            for item in newvalue[1:-1].split(','):
                nv.append(float(item))
            #print "nv=", nv
            newvalue = nv
            kw[key] = nv
            if verbose: print("newvalue=", nv, " kw=", kw)
        elif key=='epdb_flag':
            if verbose: print("setting epdb_flag to", newvalue)
            kw['epdb_flag'] = 1
        elif 'flag' in key:
            if verbose: print("key=", key, ' newvalue=', newvalue)
            if newvalue in ['1','0']:
                newvalue = int(newvalue)
            if newvalue =='False':
                newvalue = False
            if newvalue =='True':
                newvalue = True
        elif local_search and 'about' in key:
            kw['about'] = newvalue
            kw['tran0'] = newvalue
        else:
            kw[key] = newvalue
            if verbose: print("set ", key, " to ", newvalue)
        if verbose: print("calling set_docking_parameters with kw=", kw)
        dm.set_docking_parameters(*(), **kw)
        if key not in parameter_list:
            #special hacks for output_pop_file,set_sw1,...
            if key=='output_pop_file':
                parameter_list.insert(parameter_list.index('set_ga'), key)
            elif key=='set_sw1':
                parameter_list.insert(parameter_list.index('ls_search_freq')+1, key)
            else:
                parameter_list.append(key)
    if about_root_atoms_only:
        lines = dm.ligand.parser.allLines
        for ix, lll in enumerate(lines):
            if lll.find("ROOT")==0:
                root_ix = ix
            if lll.find("ENDROOT")==0:
                endroot_ix = ix
                break
        last_ix = endroot_ix - (root_ix + 1) #47-(18+1)
        crds = dm.ligand.allAtoms[0:last_ix].coords
        about = (numpy.add.reduce(crds)/float(len(crds))).tolist()
        dm.dpo['about']['value'] = (round(about[0],3), round(about[1],3), round(about[2],3))
    if epdb_output:
        dm.dpo['epdb_flag']['value'] = 1
        dm.write_dpf(dpf_filename, parm_list=epdb_list4_2)
    else:
        if verbose: print("not epdb_output")
        dm.write_dpf(dpf_filename, parameter_list, pop_seed=pop_seed)

if __name__ == '__main__':
    from pathlib import Path
    # receptor = 'C:\\Users\\o_o\\Documents\\research\\postgrad\\SOAN-UH\\E.Coli_P.Falsiparum\\structures\\03_docking_tutorial\\data\\receptor_3b2x.pdb'
    # receptor = "rec_3ebh-prep.pdbqt"
    # outreceptor = 'receptor_3b2x.pdbqt'
    # ligand = 'C:\\Users\\o_o\\Documents\\research\\postgrad\\SOAN-UH\\E.Coli_P.Falsiparum\\structures\\03_docking_tutorial\\data\\lig_3b2x.mol2'
    # ligand = "lig_3ebh-prep.pdbqt"
    # outligand = 'lig_3b2x.pdbqt'
    # outputfilename = 'C:\\Users\\o_o\\Documents\\research\\postgrad\\SOAN-UH\\E.Coli_P.Falsiparum\\structures\\03_docking_tutorial\\data\\lig_3b2x.pdbqt'
    # prepare_receptor4(receptor, verbose=True, outputfilename=outreceptor)
    # prepare_ligand4(ligand,verbose=True,charges_to_add=None, outputfilename=outligand)
    # prepare_gpf4(receptor_filename=receptor, ligand_filename = ligand, verbose=1,output_gpf_filename='here.gpf')
    # prepare_gpf4(receptor_filename=outreceptor)
    

    # python.exe "C:\Program Files (x86)\MGLTools-1.5.7rc1\Lib\site-packages\AutoDockTools\Utilities24\prepare_gpf4.py" -l C:\Users\o_o\Documents\PythonProjects\02_chemioinformatics\00_programs\mscreen\mscreen\autodocktools_prepare_py3k\AutoDockTools\Utilities24\lig_3b2x.pdbqt -r C:\Users\o_o\Documents\PythonProjects\02_chemioinformatics\00_programs\mscreen\mscreen\autodocktools_prepare_py3k\AutoDockTools\Utilities24\receptor_3b2x.pdbqt -o C:\Users\o_o\Documents\PythonProjects\02_chemioinformatics\00_programs\mscreen\mscreen\screening\file.gpf
    # python.exe "C:\Program Files (x86)\MGLTools-1.5.7rc1\Lib\site-packages\AutoDockTools\Utilities24\prepare_dpf42.py" - l C: \Users\o_o\Documents\PythonProjects\02_chemioinformatics\00_programs\mscreen\mscreen\screening\lig_3b2x.pdbqt - r C: \Users\o_o\Documents\PythonProjects\02_chemioinformatics\00_programs\mscreen\mscreen\screening\receptor_3b2x.pdbqt - o C: \Users\o_o\Documents\PythonProjects\02_chemioinformatics\00_programs\mscreen\mscreen\screening\lig_3b2x_receptor_3b2x.dpf
    # prepare_dpf42(receptor_filename=receptor,ligand_filename=ligand)