# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 19:55:51 2019

@author: lizet
"""

import os
from pathlib import Path

import numpy as np


class LogFile:
    def __init__(self, ligand, outputFolder):
        """
        The LogFile object can be used storage the vina out log file

        input
        -----
        **ligand**: the name of the ligand

        **soutputFolder**: the path to where are the result of the screening

        Examples
        --------
        screening_result/5j8o_two_chains/5j89_ligand_h/5j89_ligand_h_log.txt

        outputFolder = screening_result/5j8o_two_chains/

        ligand = 5j89_ligand_h

        >>> LogFile(ligand_01, target_x_results_folder)

        """

        self.name = ligand
        self.filepath = outputFolder / ligand
        if not self.filepath.exists():
            self.filepath = outputFolder
        for i in self.filepath.iterdir():
            if 'log' in i.name:
                log_name = i.name
        self.logPath = self.filepath / log_name
        self.GetData()
#        except ValueError:
#            print("Seems that the vina screening hasn't ended")
#

    def readLog(self):
        with open(self.logPath, 'r') as log:
            lines = log.readlines()
        return lines

    def GetData(self):
        lines = self.readLog()
        # try:
#        print(self.filepath.name)
        if 'Writing output ... done.\n' not in lines:
            self.best = np.nan
            self.ave3 = np.nan
            self.ave = np.nan
            return 'log file is not completed'
        ndx = lines.index('-----+------------+----------+----------\n')
        end = lines.index('Writing output ... done.\n')
        conformers_number = end - ndx - 1
        self.data = np.zeros(conformers_number)
        for i, j in zip(range(ndx + 1, ndx + 1 + conformers_number), range(conformers_number)):
            l = lines[i]
            self.data[j] = l.split()[1]
        self.best = self.data[0]
        self.ave3 = self.data[:3].mean()
        self.ave = self.data.mean()


class Analysis:

    def __init__(self, outputFolder="out", output_nam='result.txt'):
        self.ligands = list(d for d in os.listdir(
            outputFolder) if os.path.isdir(os.path.join(outputFolder, d)))
        self.outputFolder = Path(outputFolder)
        self.totalLigands = len(self.ligands)
        self.getLigandsResults()
        self.output_name = output_nam
        self.writeResults()

    def getLigandsResults(self):
        self.results = np.zeros((len(self.ligands), 4))
        self.ligandCode = {}
        for l, j in zip(self.ligands, range(self.totalLigands)):
            # seems to be the same as for l in enumerate(self.ligand):
            # use l[0] as j and l[1] as l
            ligData = LogFile(l, self.outputFolder)
            print(l, j)
            self.ligandCode[j] = l
            self.results[j][0] = j
            self.results[j][1] = ligData.best
            self.results[j][2] = ligData.ave3
            self.results[j][3] = ligData.ave
            ndxSorted = self.results[:, 1].argsort()
            self.results = self.results[ndxSorted]

    def writeResults(self):
        with open(self.output_name, 'w+') as result:
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
            for i in self.results:
                print(i)
        #        result.write("  {0}\t{1:.1f}\t\t{2:.2f}\t\t{3:.2f}\n".format(b.ligandCode[i[0]], i[1], i[2], i[3]))
                result.write("   {0:14s}{1:6.1f}{2:17.2f}{3:16.2f}\n".format(
                    self.ligandCode[i[0]].ljust(4), i[1], i[2], i[3]))


if __name__ == '__main__':
    print('hello world')
    # a = LogFile('JM25-SV-G02-M00', p)
