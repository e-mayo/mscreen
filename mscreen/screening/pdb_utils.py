import re
# import numpy as np
from math import sqrt

aa3 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',\
'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M',\
'ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T',\
'VAL':'V','TRP':'W','TYR':'Y','MSE':'M','GLD':'X','AuX':'X', 'UNK': '.'}



class Atom:

    def __init__(self, record, atom_id, atom_name,
                        alt_loc, res_name, chain_id,
                        res_seq, i_code, coord_x,
                        coord_y, coord_z, occ,
                        t_factor, elem, charge, string=None):
        self.record = record.strip()
        self.atom_id = int(atom_id)
        self.atom_name = atom_name.strip()
        self.alt_loc = alt_loc.strip()
        self.res_name = res_name.strip()
        self.chain_id = chain_id.strip()
        self.res_seq = res_seq.strip()
        self.i_code = i_code.strip()
        self.coord_x = float(coord_x)
        self.coord_y = float(coord_y)
        self.coord_z = float(coord_z)
        self.occ = float(occ)
        self.t_factor = float(t_factor)
        self.elem = elem.strip()
        self.charge = charge.strip()
        self.raw_string = string    
    
    def __str__(self):
        atom_str  = f'{self.record:<6}'
        atom_str += f'{self.atom_id:>5}'
        atom_str += ' '
        atom_str += f'{self.atom_name:^4}'
        atom_str += f'{self.alt_loc:^1}'
        atom_str += f'{self.res_name:>3}'
        atom_str += ' '
        atom_str += f'{self.chain_id}'
        atom_str += f'{self.res_seq:>4}'
        atom_str += f'{self.i_code}'
        atom_str += ' '*4
        atom_str += f'{self.coord_x:8.3f}'
        atom_str += f'{self.coord_y:8.3f}'
        atom_str += f'{self.coord_z:8.3f}'
        atom_str += f'{self.occ:6.2f}'
        atom_str += f'{self.t_factor:6.2f}'
        atom_str += ' '*10
        atom_str += f'{self.elem:>2}'
        atom_str += f'{self.charge:>2}'
        return atom_str

    @classmethod
    def from_str(cls,string):
        pdb_fmt = re.compile("""^(?P<record>.{6})(?P<atom_id>.{5})\s(?P<atom_name>.{4})(?P<alt_loc>.{1})(?P<res_name>.{3})\s(?P<chain_id>.{1})(?P<res_seq>.{4})(?P<i_code>.{1})\s{3}(?P<coord_x>.{8})(?P<coord_y>.{8})(?P<coord_z>.{8})(?P<occ>.{6})(?P<t_factor>.{6})\s{10}(?P<elem>.{2})(?P<charge>.{2})""",re.X|re.M)
        m = pdb_fmt.match(string)
        atom_dict = m.groupdict()
        return cls(**atom_dict)

    def to_string(self):
        return self.__str__()

    @staticmethod
    def parse_atom_from_pdb(string):
        pdb_fmt = re.compile("""^(?P<record>.{6})(?P<atom_id>.{5})\s(?P<atom_name>.{4})(?P<alt_loc>.{1})(?P<res_name>.{3})\s(?P<chain_id>.{1})(?P<res_seq>.{4})(?P<i_code>.{1})\s{3}(?P<coord_x>.{8})(?P<coord_y>.{8})(?P<coord_z>.{8})(?P<occ>.{6})(?P<t_factor>.{6})\s{10}(?P<elem>.{2})(?P<charge>.{2})""",re.X|re.M)
        m = pdb_fmt.match(string)
        return m.groupdict()



class Residue:

    def __init__(self,res_name,chain_id,res_seq,atoms=None):
        self.res_name = res_name
        self.chain_id = chain_id
        self.res_seq = res_seq
        if not atoms:
            self.atoms = []
        else:
            self.atoms = atoms

    def add_atom(self,atom):
        self.atoms.append(atom)

    def __str__(self):
        res_str  = f'{self.res_name}'
        res_str += f' {self.chain_id}'
        res_str += f' {self.res_seq}'
        return res_str
    # def get_res_as_dict(self):
    #     res = {'res_name':self.res_name,
    #            'chain_id':self.chain_id,
    #            'res_seq':self.res_seq,
    #            'atoms':self.atoms}
    #     return res



class Mol:
    pass

              


class Chain:

    def __init__(self, chain_id,residues=None, seq=None):
        self.chain_id = chain_id
        self.seq = seq
        if not residues:
            self.residues = []
        else:
            self.residues = residues

    def add_residue(self,residue):
        self.residues.append(residue)
# ^.{6}.{5}\s.{4}.{1}.{3}\s.{1}.{4}.{1}.\s{3}.{8}.{8}.{8}.{6}.{6}\s{9}.{2}.{2}

def distance(x1,y1,z1,x2,y2,z2):
    x_sqrt = (x2 -  x1)**2
    y_sqrt = (y2 -  y1)**2
    z_sqrt = (z2 -  z1)**2
    return sqrt(x_sqrt + y_sqrt + z_sqrt)

def process_pdb_file(pdb_file):
    pdb_text = read_pdb_file(pdb_file)
    pdb_fmt = re.compile("""^(?P<record>.{6})(?P<atom_id>.{5})\s(?P<atom_name>.{4})(?P<alt_loc>.{1})(?P<res_name>.{3})\s(?P<chain_id>.{1})(?P<res_seq>.{4})(?P<i_code>.{1})\s{3}(?P<coord_x>.{8})(?P<coord_y>.{8})(?P<coord_z>.{8})(?P<occ>.{6})(?P<t_factor>.{6})\s{10}(?P<elem>.{2})(?P<charge>.{2})""",re.X|re.M)
    f = pdb_fmt.finditer(pdb_text)
    res_seq = ''
    chain_id = ''
    chains = []
    res_list = []
    atoms = []
    for m in f:
        atom = Atom(**m.groupdict())
        atoms.append(atom)
        if atom.res_seq != res_seq:
            res_seq = atom.res_seq
            res = Residue(atom.res_name, atom.chain_id ,atom.res_seq)
            if atom.chain_id != chain_id:
                chain_id = atom.chain_id
                chain = Chain(chain_id)
                chains.append(chain)
            chain.add_residue(res)
            res_list.append(res)
        res.add_atom(atom)
    return chains

def get_sequence(chain):
    seq = ''
    for res in chain.residues:
        seq += aa3.get(res.res_name,'')
    return seq

def write_pdb(pdb_file,entries):
    with open(pdb_file,'w') as f:
        for entry in entries:
            f.write(f'{entry}\n')

def select(reference, radii, protein, byres=False):
    ref = [10.15, 10.08, 4.02 ] # THR`522
    ref = reference
    keep_a = []
    keep_r = set()
    # keep_a = [atom for res in protein.residues for atom in res.atoms if distance(ref[0], ref[1], ref[2], atom.coord_x, atom.coord_y, atom.coord_z) <= radii]
    for res in protein.residues:
        for atom in res.atoms:
            if distance(ref[0], ref[1], ref[2], atom.coord_x, atom.coord_y, atom.coord_z) <= radii:
                keep_a.append(atom)
                keep_r.add(atom.res_seq)   
    if byres:
        keep_r = [res for res in protein.residues if res.res_seq in keep_r]
        return [atom for res in keep_r for atom in res.atoms]
    else:
        return keep_a

# def select(reference, radii, protein,byres=False):
#     ref = [10.15, 10.08, 4.02 ] # THR`522
#     ref = reference
#     keep_a = []
#     keep_r = []
    
#     for res in protein.residues:
#         for atom in res.atoms :
#             if distance(ref[0], ref[1], ref[2], atom.coord_x, atom.coord_y, atom.coord_z) <= radii:
#                 keep_a.append(atom)
#             if byres:
#                 continue
#             keep_r.append(atom.res_seq)
#     if byres:
#         return keep_r
#     else:
#         return keep_a

def read_pdb_file(pdb_file,as_lines = False):
    with open(pdb_file) as f:
        if as_lines:
            pdb_text = f.readlines()
        else:
            pdb_text = f.read()
    return pdb_text

def extract_near_residues_as_pdb(in_pdb_file,out_pdb_file,ref,radii=20):
    """
    Read a PDB file and parse protein information.
    Select residues within radii of ref coordinates.
    Write a PDB file of the fragment after renumbering residues and atoms.
    This might fix any dms problem related to residues/atom numbers.
    Parameters:
    ----------
    in_pdb_file str:
        PDB in file name
    out_pdb_file str:
        PDB out file name
    ref list:
        Center of reference, [x, y, z]
    radii int:
        Cutoff radio. Default 20
    """
    chains = process_pdb_file(in_pdb_file)
    protein = chains[0]    
    selection = select(ref, radii, protein, byres=True)
    renumberin_atoms_list(selection)
    renumbering_residues(protein)
    write_pdb(out_pdb_file, selection)

def renumbering_residues(protein):
    for i, res in enumerate(protein.residues):
        res.res_seq = i
        for atom in res.atoms:
            atom.res_seq = i
    return protein

def renumbering_atoms(protein,n_ids=None):
    atoms = [atom for res in protein.residues for atom in res.atoms]
    for i, atom in enumerate(atoms):
        if n_ids:
            atom.atom_id =n_ids[i]
        else:
            atom.atom_id = i
    return protein

def renumberin_atoms_list(atom_list):
    for i, atom in enumerate(atom_list):
        atom.atom_id = i
    return atom_list
    
def renumbering(protein,resi=True,atomi=False, save_as=False):
    if resi: protein = renumbering_residues(protein)
    if atomi: protein = renumbering_atoms(protein)
    if save_as:
        atoms = [atom for res in protein.residues for atom in res.atoms]
        write_pdb(save_as,atoms)
    return protein
            
    

        
    
# line[0:6]  #Record name   "ATOM  "
# line[6:11]  #Integer       serial       Atom  serial number.
# line[12:16]  #Atom          name         Atom name.
# line[16]  #Character     altLoc       Alternate location indicator.
# line[17:20]  #Residue name  resName      Residue name.
# line[21]  #Character     chainID      Chain identifier.
# line[22:26]  #Integer       resSeq       Residue sequence number.
# line[26]  #AChar         iCode        Code for insertion of residues.
# line[30:38]  #Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
# line[38:46]  #Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
# line[46:54]  #Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
# line[54:60]  #Real(6.2)     occupancy    Occupancy.
# line[60:66]  #Real(6.2)     tempFactor   Temperature  factor.
# line[76:78]  #LString(2)    element      Element symbol, right-justified.
# line[78:80]  #LString(2)    charge       Charge  on the atom.
# ammino acids 3 letters code
#%%
import numpy as np


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
#%%
if __name__ == '__main__':
    
    pdb_file = 'rec_3ebh.pdb'
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
    extract_near_residues_as_pdb('rec_3ebh.pdb','out_pdb_file.pdb',[10.15, 10.08, 4.02 ],radii=20)
    