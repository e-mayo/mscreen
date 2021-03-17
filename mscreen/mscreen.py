#!/home/edd/miniconda3/envs/ai_lab/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 19:55:51 2019

@author: edd
"""
# from analysis import analysis
from screening.screen import vsprotocol
from Analysis.analysis import vsanalyser

# from analysis import analysis
import argparse
import os
import sys


def pareser():
    
    my_parser = argparse.ArgumentParser(prog="mscreen", description="Vina Screening")
    
    
    subparsers = my_parser.add_subparsers(title='task',help='The task to do. Choose between [prepare, screen, analysis]',dest='task',required=True)
    
    # virtual_screening = subparsers.add_subparsers(title='Virtual Screening',help='Virtual screening utilities')
    # my_parser2 = argparse.ArgumentParser()
    # virtual_screening = my_parser2.add_parser(title='Virtual Screening',help='Virtual screening utilities')
    
    
    
    prepare_parser = subparsers.add_parser('prepare',
                                      help='Prepare ligand and receptor for virtual screening')
    

    # subparsers.add_parser
    
    
    # prepare_parser.add_argument("-l",
    #                             "--ligands",
    #                             metavar="Raw ligands folder",
    #                             type=str,
    #                             help="The path to unprocessed ligands",
    #                             default="ligands",
    #                             required=False,
    #                             )
    
    # prepare_parser.add_argument("-r",
    #                         "--receptors",
    #                         metavar="Raw receptors folder",
    #                         type=str,
    #                         help="The path to unprocessed receptors",
    #                         default="receptors",
    #                         required=False)


    #%%
    
        parser = argparse.ArgumentParser()
        
        # subparsers = parser.add_subparsers(help='commands')
        subparsers = parser.add_subparsers(title='task',help='The task to do. Choose between [prepare, screen, analysis]',dest='task',required=True)
        # A list command
        
        list_parser = subparsers.add_parser('list', help='List contents')
        prepare_parser = subparsers.add_parser('prepare', help='Prepare ligand and receptor for virtual screening')
        prepare_parser.add_argument("-l",
                                    "--ligands",
                                    metavar="ligands_folder",
                                    type=str,
                                    help="prepare -l folder: The path to unprocessed ligands\n\
                                          screen  -l folder: The path to (un)processed ligands\
                                                             If ligands are not prepared use along with -p.",
                                    # default="ligands",
                                    required=True,
                                    )
            
        prepare_parser.add_argument("-r",
                                "--receptors",
                                metavar="receptors_folder",
                                type=str,
                                help="""prepare -r folder: The path to unprocessed receptors\n
                                      screen  -r folder: The path to (un)processed receptors\n
                                                         If receptors are not prepared use along with -p.""",
                                # default="receptors",
                                required=True)
                                # )
            
        prepare_parser.add_argument("-o",
                                "--out",
                                metavar="output_folder",
                                type=str,
                                help="Path to the output folder",
                                default="out",
                                required=False)
        
        
        
        # A create command
        screen_parser = subparsers.add_parser('screen', help='Run a virtual screening')
        
        create_parser = subparsers.add_parser('create', help='Create a directory')
        create_parser.add_argument('dirname', action='store', help='New directory to create')
        create_parser.add_argument('--read-only', default=False, action='store_true',
                                   help='Set permissions to prevent writing to the directory',
                                   )
        screen_parser.add_argument("-c",
                                    "--conf",
                                    metavar="conf_file",
                                    type=str,
                                    help="The name of the configuration file for the docking program",
                                    default="conf.txt",
                                    required=False)
         
        screen_parser.add_argument("-log",
                                "--log_file",
                                metavar="log_file",
                                type=str,
                                help="Log file name. This file anotated the time used for each docking task and some not usefull information",
                                default=None,
                                required=False)
        
        screen_parser.add_argument("-p",
                                "--prepare",
                                # metavar="autoprepare ligand and receptor",
                                # type=str,
                                help="If you will do a virtual screening with unprocessed structures use this for preprocessing",
                                action='store_true',
                                default=None,
                                required=False)
        
        analysis_parser = subparsers.add_parser('analysis', 
                                             help='Analsysis of the virtual screening output')
            
        analysis_parser.add_argument("-i",
                                    "--input",
                                    metavar="analysis_input",
                                    type=str,
                                    help="The folder containing the virtual screening result",
                                    required=True)
        
        analysis_parser.add_argument("-o",
                                    "--output",
                                    metavar="analysis_output",
                                    type=str,
                                    help="The folder for the virtual screening report",
                                    default=None,
                                    required=False)
           
        
        analysis_parser.add_argument("-t",
                            "--type_of_analysis",
                            metavar="type_of_analysis",
                            type=str,
                            help="""The type of analysis to do. 
                                    full  - return a sdf file with ligand propierties 
                                    short - write a txt file with the ligand ranking""",
                            choices=["full", "short"],
                            # usage = '%(prog)s [-a full]  -o out\
                            # return a sdf file with all the information of the \
                            #     virtual screeing\n\
                            # %(prog)s [-a short]  -o out return a txt file with\
                            #     a vina score resume',
                            default="short",
                            required=False)    
        
        parser.add_argument(
                "-d",
                "--backend",
                metavar="docking_engine",
                choices=["fwavina", "gwovina", "ledock", "psovina", "qvina2", "qvina-w", "smina", "vina", "plants"],
                type=str,
                help="The docking software to use, choose one of the next 'fwavina', 'gwovina', 'ledock', 'psovina2', 'qvina2.1', 'qvina-w', 'smina', 'vina', 'plants ",
                default="qvina-w",
                required=False,)
            
           
        parser.add_argument("-e",
                                "--exe",
                                metavar="executable_file",
                                type=str,
                                help="File containing executables path\n\
                                    eg: vina=C:\\Users\\vina.exe",
                                default=None,
                                required=False)
            
        parser.add_argument("-v", "--verbose", action="store_true")
    # my_group.add_argument('-s', '--silent', action='store_true'
#%%
    
    
    
    # prepare_parser.add_argument("-l",
    #                             "--ligands",
    #                             metavar="Ligands folder",
    #                             type=str,
    #                             help="The path to processed ligands.\
    #                                 If ligands are not prepared use along with -p.",
    #                             default="ligands",
    #                             required=False,
    #                             )
    
    # prepare_parser.add_argument("-r",
    #                         "--receptors",
    #                         metavar="Receptors folder",
    #                         type=str,
    #                         help="The path to processed receptors.\
    #                             If receptors are not prepared use along with -p.",
    #                         default="receptors",
    #                         required=False)
    
    
    
    
        
        
    # )
    # args = my_parser.parse_args()
    my_parser.add_argument("-l",
                            "--ligands",
                            metavar="ligands_folder",
                            type=str,
                            help="prepare -l folder: The path to unprocessed ligands\n\
                                  screen  -l folder: The path to (un)processed ligands\
                                                     If ligands are not prepared use along with -p.",
                            # default="ligands",
                            required=True,
                            )
    
    my_parser.add_argument("-r",
                            "--receptors",
                            metavar="receptors_folder",
                            type=str,
                            help="""prepare -r folder: The path to unprocessed receptors\n
                                  screen  -r folder: The path to (un)processed receptors\n
                                                     If receptors are not prepared use along with -p.""",
                            # default="receptors",
                            required=True)
                            # )
        
    my_parser.add_argument("-o",
                            "--out",
                            metavar="output_folder",
                            type=str,
                            help="Path to the output folder",
                            default="out",
                            required=False)

    return my_parser

parser = pareser()
parser.print_help()
#%%
if __name__ == "__main__":
    # p = pareser()
    # p
    # p.parse_args(['screen'])
    # from time import time
    # import random
    # # print(args.pkg)
    
    parser = pareser()
    args = parser.parse_args()
    #%%
    # karg = {'-l':'..\data\ligands',
    #         '-r': '..\data\receptor',
    #         '-o': '..\data\prepare',
    #         '-d': 'vina'}
    # args = parser.parse_args(['prepare',
    #                           '-l ..\data\ligands',
    #                           '-r ..\data\receptor',
    #                           '-o ..\data\prepare',
    #                           '-d vina'])
    # %%
    if not args.task:
        print("ArgumenError: task argument is required. Available task's value {prepare|screen|analysis}")
        print('Usage: python mscreen.py {prepare|screen|analysis} -l ligands_folder -r receptors_folder [-o output_folder]')
        
    print('prepare started')
    if args.task == "prepare":
        
        screening = vsprotocol.get_program(args.backend)
        s = screening(
            ligands=args.ligands,
            receptors=args.receptors,
            out=args.out,
            conf="conf.txt",#args.conf,
            backend=args.backend,
            prepare=True,
            verbose=args.verbose,
            log_file='',#args.log_file,
            exe_file=args.exe,
        )
        s.prepare_screening()

    if args.task == "screen":
        screening = vsprotocol.get_program(args.backend)
        s = screening(
            args.ligands,
            args.receptors,
            args.out,
            args.conf,
            args.backend,
            args.prepare,
            args.verbose,
            args.log_file,
            args.exe,
        )
        s.run_screening()

    elif args.task == 'analysis':
        analyser = vsanalyser.get_analyzer(args.backend)
        a = analyser(args.vs_result, args.vs_result, mode=args.type_of_analysis)
        if args.type_of_analysis == 'full':
            a.run_full_analysis()
        if args.type_of_analysis  == 'short':
            a.run_short_analysis()

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

