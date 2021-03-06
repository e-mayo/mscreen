# C:\Users\lizet\Anaconda3\python.exe
# /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 19:55:51 2019

@author: lizet
"""
# Import the argparse library

import argparse
import os
import sys
import VinaAnalysis

# Create the parser
my_parser = argparse.ArgumentParser(description='List the content of a folder')

# Add the arguments
my_parser.add_argument('-r', '--results',
                       metavar='vina results',
                       type=str,
                       help='the path to vina results',
                       default='./',
                       required=False)

my_parser.add_argument('-o', '--out',
                       metavar='out file name',
                       type=str,
                       help='the name of the out put file',
                       default='result.txt',
                       required=False)


args = my_parser.parse_args()

input_path = args.results
output_name = args.out
VinaAnalysis.Analysis(input_path, output_name)

# if __name__ == '__main__':

# input_path = '.\\test\\6r3k'
# output_name = '.\\test\\test_result_63rk_out.txt'
# VinaAnalysis.Analysis(input_path, output_name)
