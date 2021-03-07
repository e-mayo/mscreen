
#!/home/edd/miniconda3/envs/ai_lab/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 19:55:51 2019

@author: edd
"""
# from analysis import analysis
from screening.screen import vsprotocol
# from analysis import analysis
import argparse
import os
import sys


def pareser():
    my_parser = argparse.ArgumentParser(prog='mscreen',
                                        usage='%(prog)s [options]',
                                        description='Vina Screening')

    # my_group = my_parser.add_mutually_exclusive_group(required=True)

    # my_group.add_argument('-v', '--verbose', action='store_true')
    # my_group.add_argument('-s', '--silent', action='store_true')

    my_parser.add_argument('pkg',
                           metavar='pkg',
                           type=str,
                           help='screening, analysis, prepare',
                           action='store',
                           choices=['screen', 'analysis', 'prepare'])
    # usage='%(prog)s screening -l ligand_forlder -r \
    #     receptor_folder -o out  -c conf.txt /n\
    #        %(prog)s analysis  -o out -r   ')

    my_parser.add_argument('-l', '--ligands',
                           metavar='ligands folder',
                           type=str,
                           help='the path to ligands',
                           default='ligands',
                           required=False)

    my_parser.add_argument('-r', '--receptors',
                           metavar='receptors folder',
                           type=str,
                           help='the path to receptors',
                           default=None,
                           required=False)

    my_parser.add_argument('-o', '--out',
                           metavar='output_folder',
                           type=str,
                           help='The path to the output folder',
                           default='out',
                           required=False)

    my_parser.add_argument('-c', '--conf',
                           metavar='conf file',
                           type=str,
                           help='The name of the configuration file for the\
                                vina docking',
                           default='conf.txt',
                           required=False)

    my_parser.add_argument('-bk', '--backend',
                           metavar='docking_engine',
                           choices=['fwavina', 'gwovina', 'ledock', 'psovina',
                                    'qvina2', 'qvina-w', 'smina', 'vina', 'plants'],
                           type=str,
                           help="The docking software to use, choose one of the next 'fwavina', 'gwovina', 'ledock', 'psovina2', 'qvina2.1', 'qvina-w', 'smina', 'vina', 'plants ",
                           default='qvina-w',
                           required=False)

    my_parser.add_argument('-p', '--prepare',
                           metavar='ligands_folder',
                           type=str,
                           help='Prepare ligand',
                           # action='store_true',
                           default=None,
                           required=False)

    my_parser.add_argument('-rn', '--report_name',
                           metavar='report file name',
                           type=str,
                           help='Report file name',
                           default='results.txt',
                           required=False)

    my_parser.add_argument('-a', '--analysis',
                           metavar='analysis',
                           type=str,
                           choices=['full', 'short'],
                           help='The type of analysis full',
                           # usage = '%(prog)s [-a full]  -o out\
                           # return a sdf file with all the information of the \
                           #     virtual screeing\n\
                           # %(prog)s [-a short]  -o out return a txt file with\
                           #     a vina score resume',
                           default='short',
                           required=False)

    my_parser.add_argument('-log', '--log_file',
                           metavar='log_file',
                           type=str,
                           help='Log file name for further information',
                           # usage = '%(prog)s [-a full]  -o out\
                           # return a sdf file with all the information of the \
                           #     virtual screeing\n\
                           # %(prog)s [-a short]  -o out return a txt file with\
                           #     a vina score resume',
                           default=None,
                           required=False)

    my_parser.add_argument('-e', '--exe',
                           metavar='executable_file',
                           type=str,
                           help='File containing executables path',
                           # usage = '%(prog)s [-a full]  -o out\
                           # return a sdf file with all the information of the \
                           #     virtual screeing\n\
                           # %(prog)s [-a short]  -o out return a txt file with\
                           #     a vina score resume',
                           default='exe.paths',
                           required=False)

    # my_group = my_parser.add_mutually_exclusive_group(required=True)

    my_parser.add_argument('-v', '--verbose', action='store_true')
    # my_group.add_argument('-s', '--silent', action='store_true')
    args = my_parser.parse_args()

    return args


if __name__ == '__main__':

    from time import time
    import random

    args = pareser()
    if args.pkg == 'prepare':

        screening = vsprotocol.get_program(args.backend)
        s = screening(args.ligands,
                      args.receptors,
                      args.out,
                      args.conf,
                      args.backend,
                      args.prepare,
                      # args.prepared_receptors_folder,
                      # args.prepared_ligand_folder,
                      args.verbose,
                      args.log_file,
                      args.exe)
        s.prepare_screening()

    if args.pkg == 'screen':
        screening = vsprotocol.get_program(args.backend)
        s = screening(args.ligands,
                      args.receptors,
                      args.out,
                      args.conf,
                      args.backend,
                      args.prepare,
                      # args.prepared_receptors_folder,
                      # args.prepared_ligand_folder,
                      args.verbose,
                      args.log_file,
                      args.exe)
        s.run_screening()

    # elif args.pkg == 'analysis':
    #     if args.analysis == 'full':
    #         analysis.run_full_analysis(args.out)
    #     if args.analysis  == 'short':
    #         analysis.run_short_analysis(args.out)

# %%
    # class Test_mscreen:

    #     def test_prepare(self,backend=None):
    #         ligands = '../data/ePepN/ligands'
    #         receptors = '../data/ePepN/receptor'
    #         out_prepare = '../data/prepare'
    #         backend = backend
    #         screening = vsprotocol.get_program(backend)
    #         s = screening(ligands,receptors,out_prepare,backend=backend,verbose=1)
    #         s.prepare_screening()

    #     def test_screen(self,backend,config_file):
    #         ligands = f'../data/prepare/prepared_ligands_{backend}'
    #         receptors = f'../data/prepare/prepared_receptors_{backend}'
    #         out_screening = f'../data/out-{backend}'
    #         backend = backend
    #         conf = config_file
    #         screening = vsprotocol.get_program(backend)
    #         s = screening(ligands,receptors,out_screening,conf,backend=backend)
    #         s.run_screening()
# %%
    # t = Test_mscreen()
    # # t.test_prepare('vina')
    # # t.test_screen('vina','../data/config_vina_ex1.txt')
    # print('#'*144)
    # print('VINA OK')
    # print('#'*144)
    # t.test_prepare('plants')
    # t.test_screen('plants','../data/config_plants_speed4.txt')
    # print('#'*144)
    # print('PLANTS OK')
    # print('#'*144)
