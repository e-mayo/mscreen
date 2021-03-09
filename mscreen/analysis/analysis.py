# -*- coding: utf-8 -*-
"""
Created on Sun Mar 7 17:35:34 2021

@author: emayo
"""

import os
from pathlib import Path

import analysis.LogFile as LogFile
import analysis.engine as engine
import numpy as np
from rdkit import Chem


class Analysis:
    def __init__(self, input_path, output_path, mode='short'):
        """
        Analysis class

        Parameters
        ----------
        input_path : str
            The virtual screening output folder.
        output_path : str
            The output path for the analysis result.
        mode : str
            Type of analysis, short or full.

        Returns
        -------
        None.

        """
        self.input_path = Path(input_path)
        self.output_path = Path(output_path)
        self.mode = mode

    def run_short_analysis(self, filename="result.txt"):
        screening_path = Path(self.input_path)
        for rec in screening_path.iterdir():
            if not rec.is_dir():
                continue
            print(rec.name)
            self.screening_result_txt(rec)

    def run_full_analysis(self, filename="result.sdf"):
        """
        Search into each receptor folder an 
        and read all ligands output convert them to sdf
        then perfom a qt clustering of the poses using the rms
        """
        screening_path = Path(self.input_path)
        for rec in screening_path.iterdir():
            if not rec.is_dir():
                continue
            print(rec.name)
            self.screeining_result_sdf(rec)

    def write_sdf(self, file_name, mol_list, write_props=None):
        """
        Write a sdf file from a mol list

        Parameters
        ----------
        file_name : str
            The name of the sdf file.
        mol_list : list
            A list of molecules.
        write_props : list, optional
            List of propierties to writedown in the sdf file 
            eg: ['vina_pose', 'vina_score', 'cluster_id', 'clust_lenght'].
            The default is None.

        Returns
        -------
        None.

        """

        file_name = str(file_name)
        if ".sdf" not in file_name:
            file_name += ".sdf"
        w = Chem.SDWriter(file_name)
        if write_props:
            w.SetProps(write_props)
        for m in mol_list:
            #        print(m.GetProp('_Name'))
            w.write(m)
        w.close()

    def screeining_result_sdf(self, receptor_path):
        raise NotImplementedError

    def screening_result_txt(self, receptor_path, output_name="result.txt"):
        raise NotImplementedError

    def get_ligand_result(self, ligands_path):
        raise NotImplementedError


class VinaAnalysis(Analysis):
    def screeining_result_sdf(self, receptor_path):
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        threashold = 2
        site = np.array([10.2, 42.7, 40.8])
        radius = 1.5
        sdf_propierties = ["file_name", "vina_pose", "vina_score", "cluster_id", "clust_lenght", "best_pose", "in_site"]

        ligands_path = [l for l in receptor_path.iterdir() if l.is_dir()]
        for ligand_path in ligands_path:
            lig_result = self.pdbqt2sdf(ligand_path)
            lig_result = engine.qt_cluster_a_mol_list(lig_result, threashold)
            lig_result = engine.get_representative_clust(lig_result)
            for mol in lig_result:
                engine.find_mol_within_site_square(mol, site, radius, 3)
                screened_molecules[receptor_path.name].append(mol)
        sdf_name = receptor_path.parent / "screeining_result_rec_{}".format(receptor_path.name)
        self.write_sdf(sdf_name, screened_molecules[receptor_path.name], write_props=sdf_propierties)

    def screening_result_txt(self, receptor_path, output_name="result.txt"):
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        ligands_path = [l for l in receptor_path.iterdir() if l.is_dir()]

        results, ligandCode = self.get_ligand_result(ligands_path)
        out_file = self.input_path / f'screeining_result_rec_{receptor_path.name}.txt'
        self.writeResults(results, ligandCode, out_file)


    def get_ligand_result(self, ligands_path):
        totalLigands = len(ligands_path)
        results = np.zeros((len(ligands_path), 4))
        ligandCode = {}
        for l, j in zip(ligands_path, range(totalLigands)):
            # seems to be the same as for l in enumerate(ligand):
            # use l[0] as j and l[1] as l
            ligData = LogFile.LogFile(l)
            print(j, l.name)
            ligandCode[j] = l.name
            results[j][0] = j
            results[j][1] = ligData.best
            results[j][2] = ligData.ave3
            results[j][3] = ligData.ave
            ndxSorted = results[:, 1].argsort()
            results = results[ndxSorted]
        return results, ligandCode

    def writeResults(self, results, ligandCode, output_name):
        with open(output_name, "w+") as result_file:
            result_file.write(
                "\
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
# =========================================================================== #\n\n "
            )

            result_file.write(" LigandName	BestScore	AverageFT	AverageAll \n")
            for i in results:
                print(i)
                #        result.write("  {0}\t{1:.1f}\t\t{2:.2f}\t\t{3:.2f}\n".format(b.ligandCode[i[0]], i[1], i[2], i[3]))
                result_file.write(
                    "   {0:14s}{1:6.1f}{2:17.2f}{3:16.2f}\n".format(ligandCode[i[0]].ljust(4), i[1], i[2], i[3])
                )
                

        
    def run_short_analysis(self, output_nam="result.txt"):
        for rec in self.input_path.iterdir():
            if rec.is_dir():
                print(rec.name)
                self.screening_result_txt(rec)

    def run_full_analysis(self):
        """
        Search into each receptor folder
        and read all ligands-out.pdbqt convert them to sdf
        then perfom a qt clustering of the poses using the rms
        """
        for rec in self.input_path.iterdir():
            if rec.is_dir():
                print(rec.name)
                self.screeining_result_sdf(rec)

    def write_sdf(self, file_name, mol_list, write_props=None):
        """
            write an sdf storing all the molecules in a list

            input:
            -----
            file_name   --- the name of the sdf file\n
            mol_list    --- a list of molecules\n
            write_props --- list of propierties to writedown in the sdf file \
            ['vina_pose', 'vina_score', 'cluster_id', 'clust_lenght']

            example:
            --------
            write_sdf('file.sdf',[mol1,mol2],['prop1','prop2'])

        """
        file_name = str(file_name)
        if ".sdf" not in file_name:
            file_name += ".sdf"
        w = Chem.SDWriter(file_name)
        if write_props:
            w.SetProps(write_props)
        for m in mol_list:
            #        print(m.GetProp('_Name'))
            w.write(m)
        w.close()

    def sdf2mol_list(self, sdf_file):
        """
        Read a sdf file converted from vina pdbqt.

        and do some stuff

        Input:
        ------
        sdf_file: name of the sdf file

        Return:
        ------
        lig_mol_list :a list of the molecules in the sdf_file storing the
        _Name, vina_pose, vina_score
        """
        lig_mol_list = []
        sdf_file = Path(sdf_file)

        for index, mol in enumerate(Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)):
            name = mol.GetProp("_Name")  # .split('.')[0]
            mol.SetProp("_Name", name)
            try:
                pose = mol.GetProp("MODEL")
            except KeyError:
                pose = index + 1
            #            print('MODEL error at {}'.format(str(sdf_file)))
            try:
                mol.SetProp("vina_pose", str(pose))
                mol.SetProp("file_name", str(sdf_file.stem))
                pdbqtREMARK = mol.GetProp("REMARK")
                mol.SetProp("vina_score", pdbqtREMARK.split()[2])

            #                        AllChem.SanitizeMol(mol)
            except KeyError:
                pass
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True,
            )
            lig_mol_list.append(mol)
        return lig_mol_list

    def pdbqt2sdf(self, ligand_path):
        """
        convert a ligand from pdbqt to sdf
        Parameters
        ----------
        ligand_path : Path
            the ligand folder path.

        Returns
        -------
        lig_mol_list : List
            list of Chem.Mol containing all ligands poses

        """
        ligand_file = ligand_path / (ligand_path.name + "-out.pdbqt")

        if not ligand_file.exists():  # fix for the pdb bugs could be removed
            ligand_file = ligand_path / (ligand_path.name + "-out.pdb")
            # fix_pdb_bug(ligand_file)
            ligand_file = ligand_path / (ligand_path.name + "-out.pdbqt")
        lig_result = []

        sdf_file = ligand_path / (ligand_path.name + "_rdkit.sdf")
        # its required a pdbqt2sdf wrapper for drop this dependecy
        command = "obabel -ipdbqt {0} -osdf -O {1} ---errorlevel 0".format(ligand_file, sdf_file)
        if not sdf_file in (ligand_path).iterdir():
            # os.popen(command).read
            os.system(command)
        # os.popen run the program externally
        # and python flow continue and may crash
        lig_mol_list = self.sdf2mol_list(sdf_file)
        return lig_mol_list


class PlantsAnalysis(Analysis):
    pass
    
    
    
#%%
if __name__=='__main__':
    a = VinaAnalysis('../../data/out-test-vina_ex1',
                     '../../data/out-test-vina_ex1-analysis','full')
    a.run_full_analysis()
    a.run_short_analysis()