# -*- coding: utf-8 -*-
"""
Created on Wed May 26 01:07:01 2021

@author: o_o
"""

import re

# molname = 'ligand'
# replacemolname  = False


def read_file(file_name):   
    with open(file_name,'r') as f:
        # lines = f.readlines()
        text = f.read()
    return text


def get_resn_resi(text):    
    tripos = text.split('@<TRIPOS>')
    atoms = []
    for record in tripos:
        if record.startswith('ATOM'):
            lines = record.split('\n')
            for line in lines[1:-1]:
                data    = line.split()
                number     = data[0]
                oriname    = data[1]
                x          = data[2]
                y          = data[3]
                z          = data[4]
                atom_type  = data[5]
                resi       = data[6]
                resn       = data[7]
                atoms.append([resi,resn])
    return atoms

def replace(text, atoms):
    # replace ligand name
    text_r = re.sub(f'{atoms[-1][0]}\s+{atoms[-1][1]}',
                    f'{atoms[0][0]} {atoms[0][1]}',
                    text)
    # delete @<TRIPOS>SUBSTRUCTURE entry
    text_r2 = re.sub('@<TRIPOS>SUBSTRUCTURE(.*\n*)*','',text_r)
    return text_r2

def write_file(file_name, text):
    with open(file_name, 'w') as f:
        f.write(text)
        
def unique_residue(in_file, out_file):
    text = read_file(in_file)
    atoms = get_resn_resi(text)
    text_r = replace(text, atoms)
    write_file(out_file, text_r)
    

if __name__ == '__main__':
    file = 'lig_3ebh.mol2'
    unique_residue(file,f'{file}-unique.mol2') 
    
    
    
    
