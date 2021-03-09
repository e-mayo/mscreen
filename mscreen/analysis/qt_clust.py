# -*- coding: utf-8 -*-
"""
Created on Sat Jun 06 02:00:33 2020

@author: eduardo
"""
import argparse

from engine import *

my_parser = argparse.ArgumentParser(
    description='This program process the result of a vina screening done with mscreen')


my_parser.add_argument('-i', '--input',
                       metavar='input',
                       type=str,
                       help='the input sdf file containin the molecules to cluster',)

my_parser.add_argument('-o', '--output',
                       metavar='the sdf file to write',
                       type=str,
                       help='sdf_propierties to write in the sdf file',
                       default='clusterd.sdf',
                       required=False)

my_parser.add_argument('-q', '--quality_threashold',
                       metavar='quality_threashold',
                       type=float,
                       help='the quality_threashold for the qt algorithm',
                       default=2,
                       required=False)

if __name__ == "__main__":
    args = my_parser.parse_args()
    sdf_in = args.input
    sdf_out = args.output
    qt = args.quality_threashold
    mol_list = read_sdf(sdf_in)
    rmsMatrix = calc_rms_matrix(mol_list)
    mol_list = qt_cluster_a_mol_list(mol_list, qt)
    write_sdf(sdf_out, mol_list, ['vina_pose',
                                  'vina_score', 'cluster_id', 'clust_lenght'])
