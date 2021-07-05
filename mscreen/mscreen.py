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
#%%
# 

        
def pareser():
    
    main_parser = argparse.ArgumentParser(prog='mscreen',
                                          description='mScreen is an intuitive interface to many docking programs for virtual screening')
    parent_parser = argparse.ArgumentParser(add_help=False)  
    subparsers = main_parser.add_subparsers(title='Task',dest='task')

    parent_parser.add_argument(
            "-d",
            "--docking_program",
            metavar="docking_program",
            choices=["fwavina", "gwovina", "ledock", "psovina", "qvina2", "qvina-w", "smina", "vina", "plants", "dock"],
            type=str,
            help="The docking software to use, choose one of the next 'fwavina', 'gwovina', 'ledock', 'psovina2', 'qvina2.1', 'qvina-w', 'smina', 'vina', 'plants', 'dock'. Default vina ",
            default="vina",
            required=False,)
        
       
    parent_parser.add_argument("-e",
                            "--exe",
                            metavar="executable_file",
                            type=str,
                            help="File containing executables path\n\
                                eg: vina=C:\\Users\\vina.exe",
                            default=None,
                            required=False)
        
    parent_parser.add_argument("-v", "--verbose", action="store_true")
    
    prepare_parser = subparsers.add_parser('prepare',
                                      help='Prepare ligand and receptor for virtual screening',
                                      usage='mscreen prepare -l ligands_folder -r receptors_folder [-o output_folder] [-d docking_program] [-e executable_file] [-v]',
                                      parents=[parent_parser])
        
    prepare_parser.add_argument("-l",
                                "--ligands",
                                metavar="ligands_folder",
                                type=str,
                                help="prepare -l folder: The path to unprocessed ligands\n\
                                      screen  -l folder: The path to (un)processed ligands\
                                                          If ligands are not prepared use along with -p.",
                                required=True,
                                )
        
    prepare_parser.add_argument("-r",
                            "--receptors",
                            metavar="receptors_folder",
                            type=str,
                            help="""prepare -r folder: The path to unprocessed receptors\n
                                  screen  -r folder: The path to (un)processed receptors\n
                                                      If receptors are not prepared use along with -p.""",
                            required=True)
                            # )
        
    prepare_parser.add_argument("-o",
                            "--out",
                            metavar="output_folder",
                            type=str,
                            help="Path to the output folder",
                            default="out",
                            required=False)
    

    screen_parser = subparsers.add_parser('screen', 
                                          help='Run virtual screenin using [-d docking_program]',
                                          usage='mscreen screen [-h] -l ligands_folder -r receptors_folder [-o output_folder] [-c conf_file] [-log log_file] [-p] [-d docking_program] [-e executable_file] [-v]',
                                          parents=[parent_parser])
    
    screen_parser.add_argument("-l",
                                "--ligands",
                                metavar="ligands_folder",
                                type=str,
                                help="prepare -l folder: The path to unprocessed ligands\n\
                                      screen  -l folder: The path to (un)processed ligands\
                                                          If ligands are not prepared use along with -p.",
                                # default="ligands",
                                required=True,
                                )
        
    screen_parser.add_argument("-r",
                            "--receptors",
                            metavar="receptors_folder",
                            type=str,
                            help="""prepare -r folder: The path to unprocessed receptors\n
                                  screen  -r folder: The path to (un)processed receptors\n
                                                      If receptors are not prepared use along with -p.""",
                            # default="receptors",
                            required=True)
                            # )
        
    screen_parser.add_argument("-o",
                            "--out",
                            metavar="output_folder",
                            type=str,
                            help="Path to the output folder",
                            default="out",
                            required=False)

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
                                          help='Analsysis of the virtual screening results',
                                          usage='mscreen analysis [-h] -i analysis_input [-o analysis_output] [-t type_of_analysis] [-d docking_program] [-e executable_file] [-v]',
                                          parents=[parent_parser])
        
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
                        default="short",
                        required=False)    
    
    
    
    return main_parser

if __name__ == "__main__":

    parser = pareser()
    args = parser.parse_args()
    
    if not args.task:
        print("ArgumenError: task argument is required. Available task's value {prepare|screen|analysis}")
        print('Usage: python mscreen.py {prepare|screen|analysis} -l ligands_folder -r receptors_folder [-o output_folder]')
    
    
    if args.task == "prepare":
        
        screening = vsprotocol.get_program(args.docking_program)
        s = screening(
            ligands=args.ligands,
            receptors=args.receptors,
            out=args.out,
            conf="conf.txt",#args.conf,
            docking_program=args.docking_program,
            prepare=True,
            verbose=args.verbose,
            log_file='',#args.log_file,
            exe_file=args.exe,
        )
        s.prepare_screening()

    if args.task == "screen":
        screening = vsprotocol.get_program(args.docking_program)
        s = screening(
            args.ligands,
            args.receptors,
            args.out,
            args.conf,
            args.docking_program,
            args.prepare,
            args.verbose,
            args.log_file,
            args.exe,
        )
        s.run_screening()

    elif args.task == 'analysis':
        analyser = vsanalyser.get_analyzer(args.docking_program)
        a = analyser(args.input, args.input, mode=args.type_of_analysis)
        if args.type_of_analysis == 'full':
            a.run_full_analysis()
        if args.type_of_analysis  == 'short':
            a.run_short_analysis()

# # %%
class TestMscreen:

    @staticmethod
    def test_prepare(ligands="../data/ligands",
                    receptors="../data/receptors",
                    out="../data/prepare",
                    docking_program="dock" ,
                    verbose=True):
        screening = vsprotocol.get_program(docking_program)
        s = screening(ligands,receptors,out,docking_program=docking_program,verbose=1)
        s.prepare_screening()
    
    @staticmethod
    def test_screen_prepare_off(ligands="../data/prepare/prepared_ligands_",
                    receptors="../data/prepare/prepared_receptors_",
                    out="../data/out",
                    docking_program="dock",
                    conf="../data/config_dock_flex_sample.txt"):
        screening = vsprotocol.get_program(docking_program)
        ligands = ligands + docking_program
        receptors = receptors + docking_program
        s = screening(ligands,receptors,out,conf,docking_program=docking_program,prepare=False)
        s.run_screening()
    
    @staticmethod
    def test_screen_prepare_on(ligands="../data/ligands",
                    receptors="../data/receptors",
                    out="../data/out",
                    docking_program="dock",
                    conf="../data/config_dock_flex_sample.txt"):
        screening = vsprotocol.get_program(docking_program)
        s = screening(ligands,receptors,out,conf,docking_program=docking_program,prepare=True)
        s.run_screening()

# TestMscreen.test_prepare(docking_program='dock') 
# TestMscreen.test_screen_prepare_on(docking_program='dock')
# TestMscreen.test_screen_prepare_off(docking_program='dock')

