# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 07:39:31 2020

@author: lizet
"""

from numpy import array
import os

def rename_old( folder):    

    # olds_files = array([file_name for file_name in name.parent.iterdir() if name.stem in file_name])
    olds_files = []
    for i in folder.parent.iterdir():
        if folder.stem in i.name:
            olds_files.append(i)
    olds_files = array(olds_files)

    modificatios_date = array([os.path.getmtime(file_name)
                               for file_name in olds_files])
    olds_files = olds_files[modificatios_date.argsort()][::-1]
    name = folder.name
    dump_name = []
    number = 1
    for file_name in olds_files:
        file_name.rename( file_name.parent / f'd-{name}.{number}')
        # os.rename(file_name, 'd-{0}.{1}'.format(name, number))
        dump_name.append(file_name.parent / f'd-{name}.{number}')
        number += 1
    number = 1
    for file_name in dump_name:
        file_name.rename(file_name.parent / f'#{name}.{number}')
        number += 1