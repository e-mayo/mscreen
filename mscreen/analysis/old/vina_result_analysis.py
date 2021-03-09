# -*- coding: utf-8 -*-
"""
Created on Sat Jun 06 02:00:33 2020

@author: eduardo
"""

from pathlib import Path

import numpy as np
from engine import *


#%%
class VinaAnalysis():

    def run_analysis(self):
        """
        Search into each receptor folder
        and read all ligands-out.pdbqt convert them to sdf
        then perfom a qt clustering of the poses using the rms
        """
        for rec in self.screening_path.iterdir():
            if rec.is_dir():
                print(rec.name)
                self.molecules_modification(rec.name)



class VinaResults:
    """
    This class attempt to store all the information of a virutal screening
    done by using mscreen
    """

    def __init__(self, screening_out_folder, threashold=2, site=None, radius = 5, sdf_propierties=None):
        """
        initialize the VinaResults objects
        """
        if not sdf_propierties:
            self.sdf_propierties = ['vina_pose',
                                    'vina_score', 'cluster_id', 'clust_lenght']
        else:
            self.sdf_propierties = sdf_propierties
        if not site:
            self.site = np.array([10.2,42.7,40.8]) # this is just for pdl1
        else:
            self.site = site

        self.radius = radius
        self.threashold = threashold


    def run_analysis(self):
        """
        Search into each receptor folder
        and read all ligands-out.pdbqt convert them to sdf
        then perfom a qt clustering of the poses using the rms
        """
        for rec in self.screening_path.iterdir():
            if rec.is_dir():
                print(rec.name)
                self.molecules_modification(rec.name)



if __name__ == '__main__':

    print('done')
    vr = VinaResults('test\\vina_out')
    vr.run_analysis()

    f00_benchmark_crystal_exclusives_6R3K = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//02_benchmark_crystal_exclusives_6r3k//vina_result_2'
    f01_benchmark_non_exclusives_6R3K = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//02_benchmark_non_exclusives_6r3k//vina_result'
    f11_LSOS1k = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//11_LSOS1K_6r3k//vina_result'
    f11_LSOS1k_blind = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//11_LSOS1K_6r3k//vina_result_blind_docking'
    f11_JMLO_LIB_6r3k = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//11_JMLO_LIB_6r3k//vina_result'
    f11_LSO25_LIB_6r3k = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//11_LSO25_LIB_6r3k//vina_result'

    f22__lso_evolved_library_1500_6r3k = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//22__lso_evolved_library_1500_6r3k//vina_result'
    f33_em_lib_6r3k = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//33_em_lib_6r3k//vina_result'
    f34_em_169_6r3k = 'C://Users//lizet//Documents//tesis//working_on//docking//vina//pdl1//34_em_169_6r3k//vina_result'

    vr = VinaResults('test//vina_out//')
    vr.run_analysis()
# #
#     f00_benchmark_crystal_exclusives = Path('/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/02_benchmark_crystal_exclusives_6r3k/vina_result')
#     f01_benchmark_non_exclusives = Path(
#         '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/02_benchmark_non_exclusives_6r3k/vina_result')
    # f02_benchmark_crystal_exclusives_6r3k = Path(
    #     '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/02_benchmark_crystal_exclusives_6r3k/vina_result')
    # f02_benchmark_non_exclusives_6r3k = Path(
    #     '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/02_benchmark_non_exclusives_6r3k/vina_result')
    # f22__lso_evolved_library_1500_6r3k = Path(
    #     '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/22__lso_evolved_library_1500_6r3k/LSOG_LIB-6R3K-out/vina_results')
    # f32_em_169 = Path(
    #     '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/32_em_169/em-169-out')
    # f33_em_lib_6r3k = Path(
    #     '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/33_em_lib_6r3k')
    # f34_em_169_6r3k = Path(
    #     '/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/34_em_169_6r3k/E169_LIB-6R3K-out')

    paths = [f00_benchmark_crystal_exclusives_6R3K,  # 0
             f01_benchmark_non_exclusives_6R3K,  # 1
             f11_LSOS1k_blind,  # 2
             f11_LSOS1k,  # 3
             f11_JMLO_LIB_6r3k,  # 4
             f11_LSO25_LIB_6r3k,  # 5
             f22__lso_evolved_library_1500_6r3k,  # 6
             f33_em_lib_6r3k,  # 7
             f34_em_169_6r3k]  # 8
    # test = Path('/media/edd/OS/Users/lizet/Documents/tesis/working_on/docking/vina/pdl1/02_benchmark_crystal_exclusives_6r3k/vina_result_2')
    # vr = VinaResults(test)
    # vr.run_analysis()
    from time import time
    for n, i in enumerate(paths):
        i = Path(i)
        try:
            f = open('log.txt','r+')
            f.read()
            t0 = time()
            f.write('working on {}\n'.format(i.parent.name))
            print('working on {}'.format(i.parent.name))
            vr = VinaResults(i,radius=3)
            vr.run_analysis()
            f.write('dt {}\n'.format(time()-t0))
            f.close()
        except:
            f = open('log.txt','r+')
            f.read()
            t0 = time()
            f.write('error at {}\n'.format(i.parent.name))
            print('error at {}'.format(i.parent.name))
            f.write('dt {}\n'.format(time()-t0))
            f.close()

