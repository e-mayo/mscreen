# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 17:18:00 2020

@author: emayo
"""

import numpy as np
from pathlib import Path


class LogFile:
    """
    Vina Log File wrapper
    """

    def __init__(self, filepath):
        """
        init the log file
        the mscreen scripts create the next folder tree

                                  |--LIG1--|--LIG1-out.pdbqt
                                  |        |--LIG1-out-log.txt
        screen_out -|--receptor1--|--LIG2
                    |             |--LIGn
                    |
                    |             |--LIG1
                    |--receptorn--|--LIG2
                                  |--LIGn

        input --- log_file_location = scree_out/receptor/lig/
        """
        self.filepath = Path(filepath)
        for i in self.filepath.iterdir():
            if 'log' in i.name:
                log_name = i.name
        self.logPath = self.filepath / log_name
        try:
            self.GetData()
        except UnicodeDecodeError:
            print('{} UnicodeDecodeError'.format(lig_folder.name))

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