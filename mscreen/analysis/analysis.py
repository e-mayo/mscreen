# -*- coding: utf-8 -*-
"""
Created on Sun Mar 7 17:35:34 2021

@author: emayo
"""

import os
from pathlib import Path
# os.chdir('../')

import Analysis.LogFile as LogFile
import Analysis.engine as engine
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
        
    def writeResults(self, results, ligandCode, output_name):
        with open(output_name, "w+") as result_file:
            result_file.write(
                "\
# =========================================================================== #\n\
# # Results of Virtual Screening                                            # #\n\
# =========================================================================== #\n\
# This result file only contain the the ligand score                          #\n\
# The first column is the name of the ligand                                  #\n\
# The Second column is the value of the best pose                             #\n\
# The Third column is the value of the average of the first three poses       #\n\
# The fourth column is th value of the average of all poses                   #\n\
#                                                                             #\n\
#                                                                             #\n\
# =========================================================================== #\n\n "
            )

            result_file.write(" LigandName	BestScore	AverageFT	AverageAll \n")
            if not ligandCode: 
                result_file.write('empty folder')
                return None
            for i in results:
                print(i)
                #        result.write("  {0}\t{1:.1f}\t\t{2:.2f}\t\t{3:.2f}\n".format(b.ligandCode[i[0]], i[1], i[2], i[3]))
                
                result_file.write(
                    "   {0:14s}{1:6.1f}{2:17.2f}{3:16.2f}\n".format(ligandCode[i[0]].ljust(4), i[1], i[2], i[3])
                )
                
                
    @classmethod
    def is_empty(self,path):
        path = Path(path)
        if len(list(path.iterdir())) == 0:
            return True
        else:
            return False
    
    def screeining_result_sdf(self, receptor_path):
        raise NotImplementedError

    def screening_result_txt(self, receptor_path, output_name="result.txt"):
        raise NotImplementedError

    def get_ligand_result(self, ligands_path):
        raise NotImplementedError


class VinaAnalysis(Analysis):
    
    def run_short_analysis(self, output_nam="result.txt"):
        for rec in self.input_path.iterdir():
            if not rec.is_dir():continue
            print(rec.name)
            self.screening_result_txt(rec)
                
    def screening_result_txt(self, receptor_path, output_name="result.txt"):
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        ligands_path = [l for l in receptor_path.iterdir() if l.is_dir()]

        results, ligandCode = self.get_ligand_result(ligands_path)
        out_file = self.input_path / f'screeining_result_rec_{receptor_path.name}.txt'
        self.writeResults(results, ligandCode, out_file)
    

    def run_full_analysis(self):
        """
        Search into each receptor folder
        and read all ligands-out.pdbqt convert them to sdf
        then perfom a qt clustering of the poses using the rms
        """
        for rec in self.input_path.iterdir():
            if not rec.is_dir():continue
            if self.is_empty(rec):continue
            print(rec.name)
            self.screeining_result_sdf(rec)
                
    def screeining_result_sdf(self, receptor_path):
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        threashold = 2
        # site = np.array([10.2, 42.7, 40.8])
        radius = 1.5
        # sdf_propierties = ["file_name", "vina_pose", "vina_score", "cluster_id", "clust_lenght", "best_pose", "in_site"]
        sdf_propierties = ["file_name", "vina_pose", "vina_score", "cluster_id", "clust_lenght"]

        ligands_paths = [l for l in receptor_path.iterdir() if l.is_dir()]
        for ligand_path in ligands_paths:
            if self.is_empty(ligand_path):continue
            lig_result = self.pdbqt2sdf(ligand_path)
            lig_result = engine.qt_cluster_a_mol_list(lig_result, threashold)
            # lig_result = engine.get_representative_clust(lig_result,key_energy='vina_score')
            for mol in lig_result:
                # engine.find_mol_within_site_square(mol, site, radius, 3)
                screened_molecules[receptor_path.name].append(mol)
        sdf_name = receptor_path.parent / "screeining_result_rec_{}".format(receptor_path.name)
        self.write_sdf(sdf_name, screened_molecules[receptor_path.name], write_props=sdf_propierties)

    
    
    def get_log_file(self,filepath):
        log_name = None
        for i in filepath.iterdir():
            if 'log' in i.name:
                log_name = i
        return log_name
    
    def get_ligand_result(self, ligands_path):
        totalLigands = len(ligands_path)
        results = np.zeros((len(ligands_path), 4))
        ligandCode = dict()
        for l, j in zip(ligands_path, range(totalLigands)):
            # seems to be the same as for l in enumerate(ligand):
            # use l[0] as j and l[1] as l
            # if l / 
            log_name = self.get_log_file(l)
            if not log_name: continue
            ligData = LogFile.VinaLogFile(log_name)
            print(j, l.name)
            ligandCode[j] = l.name
            results[j][0] = j
            results[j][1] = ligData.best
            results[j][2] = ligData.ave3
            results[j][3] = ligData.ave
            ndxSorted = results[:, 1].argsort()
            results = results[ndxSorted]
        return results, ligandCode       
    

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
        
    
        if not ligand_file.exists():
            return []
        
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
        
        sdf_file = Path(sdf_file)
        if not sdf_file.exists():
            return []
        
        lig_mol_list = []
        
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

    
class PlantsAnalysis(Analysis):
    
    def run_short_analysis(self, output_nam="result.txt"):

        for rec in self.input_path.iterdir():
            if not rec.is_dir():continue		
            if self.is_empty(rec):continue
            print(rec.name)
            self.screening_result_txt(rec)
    
    def screening_result_txt(self, receptor_path, output_name="result.txt"):
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        ligands_path = [l for l in receptor_path.iterdir() if l.is_dir()]

        results, ligandCode = self.get_ligand_result(ligands_path)
        out_file = self.input_path / f'screeining_result_rec_{receptor_path.name}.txt'
        self.writeResults(results, ligandCode, out_file)
    
    def get_log_file(self,filepath):
        log_name = filepath / 'features.csv'
        if not log_name.exists():
            log_name = None
        return log_name
    
    def get_ligand_result(self, ligands_path):
        totalLigands = len(ligands_path)
        results = np.zeros((len(ligands_path), 4))
        ligandCode = dict()
        for l, j in zip(ligands_path, range(totalLigands)):
            # seems to be the same as for l in enumerate(ligand):
            # use l[0] as j and l[1] as l
            log_name = self.get_log_file(l)
            if not log_name: continue
            ligData = LogFile.PlantsLogFile(log_name)
            print(j, l.name)
            ligandCode[j] = l.name
            
            try : 
                results[j][0] = j
                results[j][1] = ligData.best[0]  
                results[j][2] = ligData.ave3[0]
                results[j][3] = ligData.ave[0]
            except:
                results[j][0] = j
                results[j][1] = np.nan
                results[j][2] = np.nan
                results[j][3] = np.nan
        ndxSorted = results[:, 1].argsort()
        results = results[ndxSorted]
        return results, ligandCode
    
    def run_full_analysis(self):
        """
        Search into each receptor folder
        and read all ligands-out.pdbqt convert them to sdf
        then perfom a qt clustering of the poses using the rms
        """
        for rec in self.input_path.iterdir():
            if not rec.is_dir():continue
            if self.is_empty(rec):continue
            print(rec.name)
            self.screeining_result_sdf(rec)
        
    def screeining_result_sdf(self, receptor_path):
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        threashold = 2
        ligands_paths = [l for l in receptor_path.iterdir() if l.is_dir()]
        if not ligands_paths:
            return None
        for ligand_path in ligands_paths:
            if self.is_empty(ligand_path):continue
            lig_result = self.mol22sdf(ligand_path)
            lig_result = engine.qt_cluster_a_mol_list(lig_result, threashold)
            # lig_result = engine.get_representative_clust(lig_result,
                                                         # key_energy='TOTAL_SCORE')
            for mol in lig_result:
                # engine.find_mol_within_site_square(mol, site, radius, 3)
                screened_molecules[receptor_path.name].append(mol)
        mol_list = screened_molecules[receptor_path.name]
        try:
            sdf_propierties = list(mol_list[0].GetPropsAsDict().keys())
        except IndexError:
                pass
        # lig_mol_list
        sdf_name = receptor_path.parent / "screeining_result_rec_{}".format(receptor_path.name)
        self.write_sdf(sdf_name, screened_molecules[receptor_path.name], write_props=sdf_propierties)

    def mol22sdf(self, ligand_path):
        """
        convert a ligand from mol2 to sdf using obabel
        Parameters
        ----------
        ligand_path : Path
            the ligand folder path.

        Returns
        -------
        lig_mol_list : List
            list of Chem.Mol containing all ligands poses

        """
        ligand_file = ligand_path / ("docked_ligands.mol2")


        sdf_file = ligand_path / (ligand_path.name + "_rdkit.sdf")
        # its required a pdbqt2sdf wrapper for drop this dependecy
        command = "obabel -imol2 {0} -osdf -O {1} ---errorlevel 0".format(ligand_file, sdf_file)
        if not sdf_file in (ligand_path).iterdir():
            # os.popen(command).read
            os.system(command)
        # os.popen run the program externally
        # and python flow continue and may crash
        lig_mol_list = self.sdf2mol_list(sdf_file)
        return lig_mol_list
    
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
        log_name = self.get_log_file(sdf_file.parent)
        lig_features = LogFile.PlantsLogFile(log_name)
        if not lig_features.is_file_ok: return []
        sdf_props = dict(zip(lig_features.headers, lig_features.data.T))
        

        for index, mol in enumerate(Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)):
            name = mol.GetProp("_Name")  # .split('.')[0]
            conf = int(name.split('conf_')[-1])
            
            mol.SetProp("_Name", name)
            mol.SetProp("plants pose", f"{conf}")
            for prop in sdf_props.keys():
                mol.SetProp(f"{prop}", f"{sdf_props.get(prop)[index]}")
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

class DOCK6Analysis(Analysis):
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
            print(f"analysing results for {rec.name}")
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
            print(f"analysing results for {rec.name}")
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
        else:
            sdf_propierties = list(mol_list[0].GetPropsAsDict().keys())
            w.SetProps(sdf_propierties)
        for m in mol_list:
            #        print(m.GetProp('_Name'))
            w.write(m)
        w.close()
        
    def writeResults(self, results, ligandCode, output_name):
        with open(output_name, "w+") as result_file:
            result_file.write(
                "\
# =========================================================================== #\n\
# # Results of Virtual Screening                                            # #\n\
# =========================================================================== #\n\
# This result file only contain the the ligand score                          #\n\
# The first column is the name of the ligand                                  #\n\
# The Second column is the value of the best pose                             #\n\
# The Third column is the value of the average of the first three poses       #\n\
# The fourth column is th value of the average of all poses                   #\n\
#                                                                             #\n\
#                                                                             #\n\
# =========================================================================== #\n\n "
            )

            result_file.write(" LigandName	BestScore	AverageFT	AverageAll \n")
            if not ligandCode: 
                result_file.write('empty folder')
                return None
            for i in results:
                print(i)
                #        result.write("  {0}\t{1:.1f}\t\t{2:.2f}\t\t{3:.2f}\n".format(b.ligandCode[i[0]], i[1], i[2], i[3]))
                
                result_file.write(
                    "   {0:14s}{1:6.1f}{2:17.2f}{3:16.2f}\n".format(ligandCode[i[0]].ljust(4), i[1], i[2], i[3])
                )
                
                
    @classmethod
    def is_empty(self,path):
        path = Path(path)
        if len(list(path.iterdir())) == 0:
            return True
        else:
            return False
    
    def screeining_result_sdf(self, receptor_path):
        """
        Read all folder in a mscreen/out/receptor path. And post process it. 
        Finnally writ an annotated sdf file.
    
        Parameters
        ----------
        receptor_path : pathlib.Path
            Directory of the virtual screening output.
    
        Returns
        -------
        None.
    
        """
        
        receptor_path = Path(receptor_path)
        screened_molecules = {}
        screened_molecules[receptor_path.name] = []
        
        threashold = 2  
        sdf_propierties = ['Name',
                         'DOCK_Rotatable_Bonds',
                         'Molecular_Weight',
                         'Formal_Charge',
                         'Cluster_Size',
                         'Grid_Score',
                         'Grid_vdw_energy',
                         'Grid_es_energy',
                         'Internal_energy_repulsive']
        sdf_propierties = None
        
        
        
        ligands_paths = [l for l in receptor_path.iterdir() if l.is_dir()]
        
        for ligand_path in ligands_paths:
            if self.is_empty(ligand_path):continue
            for file in ligand_path .iterdir():
                if 'scored' in str(file):
                    ligand_scored = file 
                    print(ligand_scored)
            lig_result = self.result2sdf(ligand_scored)
            lig_result = engine.qt_cluster_a_mol_list(lig_result, threashold)
            mol_props, mol_blocks = self.read_dock6_mol2fie(ligand_scored)
            # lig_result = engine.get_representative_clust(lig_result,key_energy='vina_score')
            for mol in lig_result:
                # engine.find_mol_within_site_square(mol, site, radius, 3)
                screened_molecules[receptor_path.name].append(mol)
        sdf_name = receptor_path.parent / "screeining_result_rec_{}".format(receptor_path.name)
        self.write_sdf(sdf_name, screened_molecules[receptor_path.name], write_props=sdf_propierties)
     
    @classmethod
    def result2sdf(self,file):
        mol_props, mol_blocks = self.read_dock6_mol2fie(file)
        lig_mol_list = []
        for index, mol_block in enumerate(mol_blocks):
            # this shoud be refined becouse some molecules may fail.
            mol = Chem.MolFromMol2Block(mol_block) 
            if not mol: continue
            props = mol_props[index]
            for key, value in props.items():
                if key == "Name":
                    mol.SetProp("_Name", value)
                else:
                    mol.SetProp(key, value)
                
                mol.SetProp("dock6_pose", str(index + 1))                
            lig_mol_list.append(mol)
        return lig_mol_list
    
    @classmethod
    def read_dock6_mol2fie(self,file):    
        with open(file,'r') as f:
            text = f.read()
            # molblock_text = re.sub("#.*\n", "", text)
            # props_text = re.sub("^(?!#).+$", "", text)
            # props_text = re.sub("^(?<!#).+$", "", text,re.MULTILINE)
        mol_props = self.get_mol_props(text)
        mol_blocks = self.get_mol_blocks(text)
        return mol_props, mol_blocks
    
    @classmethod    
    def get_mol_props(self,text):
        mols_props = []
        mol_props = {}
        for line in text.split('\n'):
            if line.startswith("@<TRIPOS>MOLECULE"):
                mols_props.append(mol_props)
                mol_props = {}
            if not line.startswith('##########'): continue
            prop = line.replace('##########','').split(':')
            mol_props.setdefault(prop[0].strip(),prop[1].strip())
            
        return mols_props
    
    @classmethod
    def get_mol_blocks(self,text):
        mol_blocks = []
        mol_block = ''
        for line in text.split('\n'):
            if line.startswith('##########'): continue
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol_blocks.append(mol_block)
                mol_block = ''
            mol_block += f"{line}\n"
        mol_blocks = self._check_mol_block_list(mol_blocks)
        return mol_blocks
    
    @classmethod
    def _check_mol_block_list(self,mol_block_list):
        block_start = "@<TRIPOS>MOLECULE"
        cleaned_list = [mol_block for mol_block in mol_block_list if block_start in mol_block]
        return cleaned_list
    
    @classmethod
    def get_log_file(self,filepath):
        log_name = None
        for i in filepath.iterdir():
            if 'score' in i.name:
                log_name = i
        return log_name
    
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
        ligandCode = dict()
        for l, j in zip(ligands_path, range(totalLigands)):
            # seems to be the same as for l in enumerate(ligand):
            # use l[0] as j and l[1] as l
            # if l / 
            log_name = self.get_log_file(l)
            if not log_name: continue
            ligData = LogFile.DockLogFile(log_name)
            print(j, l.name)
            ligandCode[j] = l.name
            results[j][0] = j
            results[j][1] = ligData.best
            results[j][2] = ligData.ave3
            results[j][3] = ligData.ave
            ndxSorted = results[:, 1].argsort()
            results = results[ndxSorted]
        return results, ligandCode       


class AnalysisFactory:
    """
    Config file reader object. This is a factory-like class.
    It have 3 functions:
        register reader
        get_reader
        read_conf       
    You shoud create an ConfReader object then register each reader you want
    and then use read_conf to read the conf file. 
    """

    def __init__(self):
        self._analysers = {}

    def register_analyzer(self, docking_program, analyzer):
        """
        This function add an analyser to the _analysers dictionary

        Parameters
        ----------
        docking_program : str
            Key of the reader.
        screening : Screening
            Screening engine to be used(Should inherit from Screening class).

        Returns
        -------
        None.

        """
        self._analysers[docking_program] = analyzer

    def get_analyzer(self, docking_program):
        """
        Get screening object

        Parameters
        ----------
        docking_program : str
            Reader key.

        Raises
        ------
        ValueError
            If docking_program dont match any _readers keys.

        Returns
        -------
        reader : Reader
            Return the reader object.

        """
        analyzer = self._analysers[docking_program]
        if not analyzer:
            raise ValueError(docking_program)
        return analyzer


    def run_analysis(self, result_folder,atype, docking_program, **kwargs):
        """
        

        """
        analyzer = self.get_analyzer(docking_program)
        a = analyzer(**kwargs)
        a = VinaAnalysis(result_folder,result_folder,atype)
        if atype =='short':
            a.run_short_analysis()
        elif atype =='full':
            a.run_full_analysis()    
        else: 
            print("Analysis type must be 'short' or 'full'")
        return 1

vsanalyser = AnalysisFactory()
vsanalyser.register_analyzer('vina', VinaAnalysis)
vsanalyser.register_analyzer('plants', PlantsAnalysis)
vsanalyser.register_analyzer('dock', DOCK6Analysis)
#%%
if __name__=='__main__':    
    # a = VinaAnalysis('../data/out-test-vina_ex1',
    #                    '../data/out-test-vina_ex1-analysis','full')
    # a.run_full_analysis()
    
    # a.run_short_analysis()
    # p = PlantsAnalysis('C:/Users/o_o/Documents/research/postgrad/SOAN-UH/E.Coli_P.Falsiparum/VS/00_results_raw/ePepN-out/out-plants_speed1','C:/Users/o_o/Documents/research/postgrad/SOAN-UH/E.Coli_P.Falsiparum/VS/00_results_raw/ePepN-out/out-plants_speed1','short')
    # # p = PlantsAnalysis('../data/out-test-vina_ex1',
    #                   # '../data/out-test-plants_speed4','full')
    # p.run_short_analysis()
    # p.run_full_analysis()
    # vs_out = 'C:/Users/o_o/Documents/research/postgrad/SOAN-UH/E.Coli_P.Falsiparum/VS/00_results_raw/ePepN-out'
    vs_out = '../../data'
    vs_out = Path(vs_out)
    for folder in vs_out.iterdir():
        if not folder.is_dir():continue
        if 'ina' in folder.name:
            a = VinaAnalysis(folder,folder,'short')
            print('running short analysis')
            a.run_full_analysis()    
            print('running full analysis')
            a.run_short_analysis()
        elif  'plants' in folder.name:
            p = PlantsAnalysis(folder,folder,'short')
            print('running short analysis')
            p.run_full_analysis()    
            print('running full analysis')
            p.run_short_analysis()
        elif  'dock' in folder.name:
            p = DOCK6Analysis(folder,folder,'short')
            print('running short analysis')
            p.run_full_analysis()    
            print('running full analysis')
            p.run_short_analysis()
