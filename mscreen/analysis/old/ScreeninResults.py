# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 17:42:22 2020

@author: lizet
"""
import re
import numpy as np
import os
from pathlib import Path
from engine import *
from LogFile import LogFile

class ScreeninResults():
    """
    All the screeining results are here.
    """
    def __init__(self, screening_out_folder):
        """
        Take the screeining folder and search into it all the information

        Parameters
        ----------
        screeining_out_folder : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.screening_path = Path(screening_out_folder)
        self.receptors_path = {} # keys are the receptor name
        self.ligands = {} #
        self.results = {}
        self.ligandCode = {}

        for rec in self.screening_path.iterdir():
            if rec.is_dir():
                print(rec.name)
                self.receptors_path[rec.name] = rec
                ligs = []
                for lig in rec.iterdir():
                    if lig.is_dir():
                        ligs.append(lig.name)
                self.ligands[rec.name] = ligs
                self.get_receptor_results(rec.name)
                self.write_receptor_results(rec.name)
        self.num_receptors = len(self.receptors_path)

    def get_receptor_results(self, receptor):
        """
        Look into all the ligands inside a receptor folder
        then look for each ligand log file and take read it

        Parameters
        ----------
        receptor : str
            name of the receptor.

        """

        ligands = self.ligands[receptor]
        self.results[receptor] = np.zeros((len(ligands), 4))
        self.ligandCode[receptor] = {}
        # for l, j in zip(self.ligands, range(self.totalLigands)):
        for j, l in enumerate(self.ligands[receptor]):
            # should be use enumerate insted
            # seems to be the same as for l in enumerate(self.ligand):
            # use l[0] as j and l[1] as l
            lig_folder = self.receptors_path[receptor] / l
            try:
                ligData = LogFile(lig_folder)
            except:
                print(lig_folder.name)
                print('error at ligand {} log_file_not_found'.format(lig_folder.name))
#            print(l, j)
            self.ligandCode[receptor][j] = l
            self.results[receptor][j][0] = j
            self.results[receptor][j][1] = ligData.best
            self.results[receptor][j][2] = ligData.ave3
            self.results[receptor][j][3] = ligData.ave
            ndxSorted = self.results[receptor][:, 1].argsort()
            self.results[receptor] = self.results[receptor][ndxSorted]
        return self.ligandCode[receptor], self.results[receptor]

    def write_receptor_results(self, receptor, file_name='result.txt'):
        """
        create a filer result.txt and write the receptor result
        input: --- file_path

        """
        file = self.receptors_path[receptor] / file_name
        with open(file, 'w+') as result:
            result.write("\
# =========================================================================== #\n\
# # Results of Vina Screening                                               # #\n\
# =========================================================================== #\n\
#                                                                             #\n\
# The first column is the name of the ligand                                  #\n\
# The Second column is the value of the best pose                             #\n\
# The Third column is the value of the average of the first three poses       #\n\
# The fourth column is th value of the average of all poses                   #\n\
#                                                                             #\n\
#                                                                             #\n\
# =========================================================================== #\n\n ")

            result.write(" LigandName	BestScore	AverageFT	AverageAll \n")
            for i in self.results[receptor]:
                #                print(i)
                #        result.write("  {0}\t{1:.1f}\t\t{2:.2f}\t\t{3:.2f}\n".format(b.ligandCode[i[0]], i[1], i[2], i[3]))
                result.write("   {0:14s}{1:6.1f}{2:17.2f}{3:16.2f}\n".format(
                    self.ligandCode[receptor][i[0]].ljust(4), i[1], i[2], i[3]))


if __name__ == '__main__':

    import argparse
    import os
    import sys
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
    sr = ScreeninResults(input_path)
