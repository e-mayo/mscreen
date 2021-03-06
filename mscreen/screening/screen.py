# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 22:21:20 2021

@author: o_o
"""

import os
import sys
import re
from time import time
from pathlib import Path
from shutil import copy, rmtree, which


# from utils import rename_old
# import prepare_vina
# import prepare_dock
# # from utils import multi_threading, count_processors

# from config import reader
# from config import writer

from screening.utils import rename_old
from screening import prepare_vina
from screening import prepare_dock
from multiprocess.utils import multi_threading, count_processors

from screening.config import reader
from screening.config import writer


class ScreeningFactory:
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
        self._program = {}

    def register_program(self, docking_program, screening):
        """
        This function add a reader to the _readers dictionary

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
        self._program[docking_program] = screening

    def get_program(self, docking_program):
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
        screening = self._program[docking_program]
        if not screening:
            raise ValueError(docking_program)
        return screening

    def prepare_screening(self, file_name, docking_program, **kwargs):
        screening = self.get_program(docking_program)
        s = screening(**kwargs)
        s.prepare_screening()
        return 1

    def run_vscreening(self, file_name, docking_program, **kwargs):
        """
        If you want to read a conf file use this function.
        It read the configuration file usinng the docking_program register 
        in _readers.


        Parameters
        ----------
        file_name : Path
            The conf file path.
        docking_program : str
            The docking_program to use. eg: 'vina'

        Returns
        -------
        Dict
            A dictionary containing the information of the conf file.

        Example:
        --------
        I: readeer.read_conf('conf.txt','vina')
        O: {'receptor': 'arg',
             'flex': 'arg',
             'ligand': 'arg',
             'center_x': 'arg',
             'center_y': 'arg'}

        """
        screening = self.get_program(docking_program)
        s = screening(**kwargs)
        s.run_screening()

        return 1


class Screening:
    """
    Templete class for all program-specific screening class 
    """

    def __init__(
        self,
        ligands,
        receptors,
        out="out",
        conf="conf.txt",
        docking_program="vina",
        prepare=False,
        verbose=False,
        log_file=None,
        exe_file="exe.paths",
    ):

        self.SUPPORTED_FORMATS = (
            ".mol2", ".pdbqt", ".pdb", ".pdbqt", ".cif", ".pqr")

        self.verbose = verbose
        self.docking_program = docking_program
        self.ligand_folder = Path(ligands)
        self.ligands = list(self.ligand_folder.iterdir())
        # self.ligands = []
        self.ligands_pdbqt = []
        self.receptors_pdbqt = []
        self.prepare = prepare
        self.receptors_prepared = []
        self.ligands_prepared = []
        self.receptors_folder = Path(receptors)
        self.receptors = list(self.receptors_folder.iterdir())
        # self.receptors = []
        self.out_folder = Path(out)
        self.conf = Path(conf)
        # self.prepared_ligand_folder = prepared_ligand_folder
        # self.prepared_receptors_folder = prepared_receptors_folder
        self.log_file = log_file
        if exe_file:
            self.exe = self.read_exeFile(Path(exe_file))

    def read_exeFile(self, file):
        exe_dict = dict()
        with open(file, "r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                line = line.replace("\n", "")
                program, path = line.split("=")
                exe_dict.setdefault(program, path)
        return exe_dict

    def logger(self, line):
        if not self.log_file:
            return None
        if self.verbose:
            print(line)
        self.log_file = Path(self.log_file)
        if not self.log_file.exists():
            with open(self.log_file, "w") as f:
                f.write("#Virtual screening logfile:\n")
                f.write(f"#docking_program: {self.docking_program}\n")
                f.write(f"#conf: {self.conf}\n")
                f.write(f"#ligand_folder: {self.ligand_folder}\n")
                f.write(f"#receptors_folder: {self.receptors_folder}\n")
                f.write("program,receptor,ligand,time\n")

        with open(self.log_file, "r+") as f:
            f.read()
            f.write(f"{line}\n")
        return line

    def prepared_folder(self, folder_name=None):
        """
        This function create the folder for the prepared files if they doesn't exists
        """
        if not self.out_folder.exists():
            self._create_output_folder_()

        if folder_name:
            folder_name = "_" + folder_name
        else:
            folder_name = ""

        self.prepared_ligand_folder = self.out_folder / \
            f"prepared_ligands{folder_name}"
        if not self.prepared_ligand_folder.exists():
            self.prepared_ligand_folder.mkdir()

        self.prepared_receptors_folder = self.out_folder / \
            f"prepared_receptors{folder_name}"
        if not self.prepared_receptors_folder.exists():
            self.prepared_receptors_folder.mkdir()

    def init_screening(self):
        """
        This function create the outpout directory tree.
        If there is an output directory with the same name rename it as 
        #{output_name}.N Where n is an integer the bigger the older
        is the folder 

        Returns
        -------
        None.

        """

        self._create_output_folder_()
        self.receptors_name = []
        self.ligands_name = []
        self.ligands = []  # this is a quick fix
        self.receptors = []  # this is a quick fix

        # these for loops should be extracted as method
        for rec in self.receptors_folder.iterdir():
            if not self._check_format_(rec):
                continue
            if rec.stem in self.receptors_name:
                continue
            self.receptors.append(rec)
            self.receptors_name.append(rec.stem)

            os.mkdir(self.out_folder / rec.stem)

        for lig in self.ligand_folder.iterdir():
            if not self._check_format_(lig):
                continue
            if lig.stem in self.ligands_name:
                continue
            self.ligands.append(lig)
            self.ligands_name.append(lig.stem)

            # out_path = self.out_folder / rec.stem / lig.stem
            # os.mkdir(out_path)

    def _check_format_(self, file):
        if file.suffix in self.SUPPORTED_FORMATS:
            return True
        else:
            print(f"{file.name} has not a supported formats")
            print(
                f"Please convert it to any supported \
                  formats: {self.SUPPORTED_FORMATS}"
            )
            return False

    def _create_output_folder_(self):
        """
        Check if the outputfolder already exists if so rename it 
        and create outputfolder
        """
        if self.out_folder.exists():
            rename_old(self.out_folder)
            os.mkdir(self.out_folder)
        else:
            os.mkdir(self.out_folder)

    def prepare_screening(self):
        """
        This function prepare ligands and receptor with the required format 
        of the docking program used. It must be overwrited 
        in all screening class.         
        """
        raise NotImplementedError

    def run_screening(self):
        """
        This function run the virtual screening all results must be writed
        in the directory tree create by init_screening.It must be overwrited 
        in all screening class.         
        """
        raise NotImplementedError

    def check_program_executable(self, docking_executable):
        """
        Check if the docking program executable exists.
        """
        if docking_executable.exists() is False:
            printout = "Docking executable could not be found at: "
            printout = printout + "{}".format(docking_executable)
            print(printout)
            raise Exception(printout)

    def _get_excutable_path_(self, exe_name):
        # from distutils.spawn import find_executable
        if which(exe_name):
            return exe_name
        else:
            return False

    def _get_executable_user_(self, exe_name):
        if not self.exe:
            return None
        exe_path = Path(self.exe.get(exe_name))
        if not exe_name:
            return False
        if exe_path.exists():
            return exe_path
        else:
            return False

    def _get_executable_dist_(self, exe_name):

        bin_path = Path(os.path.realpath("docking_executable")
                        ).parent / "screening" / "docking_executable"

        linux_path = bin_path / "linux" / exe_name
        windows_path = bin_path / "win" / exe_name

        if sys.platform == "linux" or sys.platform == "linux2":
            docking_executable = linux_path
            self.check_program_executable(docking_executable)
            return docking_executable

        elif sys.platform == "win32":
            docking_executable = windows_path
            self.check_program_executable(docking_executable)
            # print(docking_executable)
            return docking_executable

        else:
            return False

    def _get_executable_(self, exe_name):
        """
        Get the rigth executable.

        Parameters
        ----------
        windows_path : pathlib.Path
            DESCRIPTION.
        linux_path : pathlib.Path
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        docking_executable : TYPE
            DESCRIPTION.

        """
        exe = self._get_excutable_path_(exe_name)
        if exe:
            if self.verbose:
                print(f"Using {exe}")
            return exe
        exe = self._get_executable_user_(exe_name)
        if exe:
            if self.verbose:
                print(f"Using {exe}")
            return exe
        exe = self._get_executable_dist_(exe_name)
        if exe:
            if self.verbose:
                print(f"Using {exe}")
            return exe
        raise FileNotFoundError("Not excecutable found")

    def fix_multiplefragment(self, mol_file):
        with open(mol_file, "r+") as f:
            mol_text = f.read()
            p = re.compile("(\d\s\s)(\w{2,})", re.M | re.I)
            # lig_code = p.search(mol_text)[2]
            lig_code = "LIG"
            nmol_text = p.sub(f"\g<1>{lig_code}", mol_text)
        with open(mol_file, "w") as f:
            f.write(nmol_text)
        return nmol_text

#

class DOCKScreening(Screening):
    def __init__(
        self,
        ligands,
        receptors,
        out="out",
        conf="conf.txt",
        docking_program="dock",
        prepare=False,
        verbose=False,
        log_file=None,
        exe_file="exe.paths",
    ):

        self.SUPPORTED_FORMATS = (
            ".mol2", ".pdb")

        self.verbose = verbose
        self.docking_program = docking_program
        self.ligand_folder = Path(ligands)
        self.ligands = list(self.ligand_folder.iterdir())
        # self.ligands = []
        self.ligands_mol2 = []
        self.receptors_mol2 = []
        self.prepare = prepare
        self.receptors_prepared = []
        self.ligands_prepared = []
        self.receptors_folder = Path(receptors)
        self.receptors = list(self.receptors_folder.iterdir())
        # self.receptors = []
        self.out_folder = Path(out)
        self.conf = Path(conf)
        # self.prepared_ligand_folder = prepared_ligand_folder
        # self.prepared_receptors_folder = prepared_receptors_folder
        self.log_file = log_file
        if exe_file:
            self.exe = self.read_exeFile(Path(exe_file))

    def prepare_receptors(self):
        """
        Prepare receptor for DOCK6. The receptor must be preprocessed 
        (delete solvent, add hiydrogens, fill missing residues/atoms). 
        Prepare means:
        - generate receptor surface  [dms]
        - calculate spheres  [sphgen] 
        - select sphere within x Angstrom from the binding site [sphere_selector]
        - box generation [showbox]
        - grid generation [grid]
        """
        
        receptors_prepared = []
        for rec in self.receptors:
            # pdb and mol2 supported by DOCK, 
            # but extract nearatoms is not implemented for pdb file
            self.SUPPORTED_FORMATS = ['.pdb']
            if not self._check_format_(rec):
                continue
            if self.verbose:
                print(f"\n\npreparing receptor {rec.name}")
            out = self.prepared_receptors_folder / f"{rec.stem}"
            
            # preprocess receptor using prepdock # only output as mol2, PDB is needed
            # prepare_dock.dockprep_chimera(rec, out)
            
            # generating required file for docking (dms, sph, grid)
            prepare_dock.prepare_receptor_dock(rec,
                                               self.ligands[0],
                                               out_folder = self.prepared_receptors_folder,
                                               sph_selector_cut_off=10.0,
                                               ref=None) 
            
            receptors_prepared.append(out)
        self.receptors = receptors_prepared

    def prepare_ligands(self):
        """
        This function prepare ligands for DOCK6
        """
        ligands_prepared = []
        self.SUPPORTED_FORMATS = ['.mol2']
        for lig in self.ligands:
            # .mol2
            if not self._check_format_(lig):
                continue
            if self.verbose:
                print(f"\n\npreparing ligand {lig.name}")
            # self.fix_multiplefragment(lig)
            kwargs = {"verbose": True}
            out = self.prepared_ligand_folder / f"{lig.stem}.mol2"
            prepare_dock.prepare_ligand_dock(lig, out)
            ligands_prepared.append(out)
        self.ligands = ligands_prepared

    def prepare_screening(self):
        """
        This function prepare ligands and receptor with the required format 
        of the docking program used. It must be overwrited 
        in all screening class.         
        """
        self.prepared_folder(folder_name="dock")
        self.prepare_receptors()
        self.prepare_ligands()

    def run_screening(self, init_screening=True):

        if init_screening:
            self.init_screening()

        if self.prepare:
            self.logger("#Preparing ligand and receptors")
            self.prepare_screening()

        self.logger("#" * 72)
        self.logger(f"#Running virtual screening using {self.docking_program}")
        self.logger("#" * 72)

        for rec in self.receptors:
            # if rec.suffix != ".mol2":
                # continue
            list_args = []
            for lig in self.ligands:
                # if lig.suffix != ".mol2":
                    # continue
                if self.verbose:
                    self.logger(f"#Docking {lig.name} in {rec.name}")
                #                if os.listdir(self.ligand_folder).index(lig) <= x:
                #                    continue
                out_path = self.out_folder / rec.stem / lig.stem
                os.mkdir(out_path)
                t0 = time()
                self.run_dock(out_path, lig, rec)
                list_args.append((out_path, lig, rec))
                self.logger(
                    f"{self.docking_program},{rec.stem},{lig.stem},{time()-t0:.2f}")
    
    def run_dock(self, out_path, lig_path, rec_path):

        conf_file = self.out_folder / rec_path.stem / \
            f"{rec_path.stem}-{lig_path.stem}-conf.txt"
        
        rec_file = rec_path
        rec_name = rec_file.stem
        
        if 'frag-' in rec_name: # if the recepetor is already prepared skip fragment renaming
            rec_name = rec_name.replace('frag-','')
        
        frag_name = f"frag-{rec_name}"

        frag_file = rec_file.parent / f"{frag_name}.pdb"
        sph_file = rec_file.parent / f"{frag_name}.sph"
        sel_file = rec_file.parent / f'{sph_file.stem}-sele.sph'
        grid_prefix = rec_file.parent / f'{rec_file.stem}-grid'.replace('frag-','') # solve frag-naming
        ligand_outfile_prefix = out_path / f'{rec_file.stem}-{lig_path.stem}-flex-out'
        
        keywords = {
            "receptor_site_file": sel_file, #sph_file
            "ligand_atom_file": lig_path,
            "grid_score_grid_prefix": grid_prefix,  #grid_prefix.
            "ligand_outfile_prefix": ligand_outfile_prefix,
        }

        # print(reader.read_conf(self.conf, "dock"))
        dock_keywords = reader.read_conf(self.conf, self.docking_program)
        dock_keywords = {**dock_keywords, **keywords}
        writer.write_conf(dock_keywords, "dock", conf_file)
        docking_executable = self.get_docking_executable(self.docking_program)
        command = f"{docking_executable} -i {conf_file} -o {conf_file}-out.txt "
        os.system(command)
        # os.popen(command).read()

    
    def get_docking_executable(self, docking_program):

        docking_engine = {
            "dock": "dock6",  # linux & winc
        }

        exe_name = docking_engine[docking_program]
        return self._get_executable_(exe_name)

class VinaScreening(Screening):
    def prepare_receptors(self):
        """
        This function prepare receptors for autodock vina
        """
        receptors_prepared = []
        for rec in self.receptors:
            # pdb,mol2,pdbq,pdbqs,pdbqt, possibly pqr,cif
            if not self._check_format_(rec):
                continue
            if self.verbose:
                print(f"\n\npreparing receptor {rec.name}")
            kwargs = {"verbose": True}
            out = self.prepared_receptors_folder / f"{rec.stem}.pdbqt"
            prepare_vina.prepare_receptor_vina(rec, out, **kwargs)
            receptors_prepared.append(out)
        self.receptors = receptors_prepared

    def prepare_ligands(self):
        """
        This function prepare ligands for autodock vina
        """
        ligands_prepared = []
        for lig in self.ligands:
            # .pdb .mol2 or .pdbqt
            if not self._check_format_(lig):
                continue
            if self.verbose:
                print(f"\n\npreparing ligand {lig.name}")
            # self.fix_multiplefragment(lig)
            kwargs = {"verbose": True}
            out = self.prepared_ligand_folder / f"{lig.stem}.pdbqt"
            prepare_vina.prepare_ligand_vina(lig, out, **kwargs)
            ligands_prepared.append(out)
        self.ligands = ligands_prepared

    def prepare_screening(self):
        """
        Look into receptors folder and ligands folder and prepare each
        file for docking. This use a modification of prepare_ligand4
        and prepare_receptor4 so they can be called as functions.

        Parameters
        ----------
        verbose : bool, optional
            Print what its doing. The default is False.

        Returns
        -------
        None.

        """
        self.SUPPORTED_FORMATS = (
            ".mol2", ".pdbqt", ".pdb", ".pdbqt", ".cif", ".pqr")
        self.prepared_folder(folder_name="vina")
        self.prepare_receptors()
        self.prepare_ligands()

    def get_docking_executable(self, docking_program):

        docking_engine = {
            "fwavina": "fwavina",  # linux & win
            "vina": "vina",  # linux & win
            "gwovina": "gwovina",  # linux
            "ledock": "ledock",  # linux
            "psovina": "psovina2",  # linux
            "qvina2": "qvina2.1",  # linux
            "qvina-w": "qvina-w",  # linux
            "smina": "smina",  # linux
        }

        exe_name = docking_engine[docking_program]
        return self._get_executable_(exe_name)

    def run_vina(self, out_path, lig_path, rec_path):

        if self.conf.exists() == False:
            print("Configuration file is required for virtual screening")
            print("edning docking process")
            return None
            # raise ValueError
        copy(self.conf, self.out_folder / "conf.txt")
        log = out_path / "{0}-log.txt".format(lig_path.stem)
        out = out_path / "{}-out.pdbqt".format(lig_path.stem)
        docking_executable = self.get_docking_executable(self.docking_program)

        command = f"{docking_executable} --config {self.conf}\
                    --receptor {rec_path}\
                    --ligand {lig_path}\
                    --out {out} --log {log}"
        os.system(command)
        # os.popen(command).read()

    def run_screening(self, init_screening=True):

        self.logger("#" * 72)
        self.logger(f"#Running virtual screening using {self.docking_program}")
        self.logger("#" * 72)

        if init_screening:
            self.init_screening()

        if self.prepare:
            self.logger("#Preparing ligand and receptors")
            self.prepare_screening()

        for rec in self.receptors:
            if rec.suffix != ".pdbqt":
                continue
            for lig in self.ligands:
                if lig.suffix != ".pdbqt":
                    continue
                self.logger(f"#Docking {lig.name} in {rec.name}")
                #                if os.listdir(self.ligand_folder).index(lig) <= x:
                #                    continue
                out_path = self.out_folder / rec.stem / lig.stem
                os.mkdir(out_path)
                t0 = time()
                self.run_vina(out_path, lig, rec)
                self.logger(
                    f"{self.docking_program},{rec.stem},{lig.stem},{time()-t0:.2f}")

    # def set_metalloprotein_charge(self,atomname,charge):
    # pass


class PlantsScreening(Screening):
    def prepare_screening(self):
        """
        Look into receptors folder and ligands folder and prepare each
        file for docking. This use a SPORE for preparing ligand and receptor

        Parameters
        ----------
        verbose : bool, optional
            Print what its doing. The default is False.

        Returns
        -------
        None.

        """
        self.SUPPORTED_FORMATS = [".mol2", ".pdb"]
        self.prepared_folder(folder_name="plants")
        self.prepare_receptors()
        self.prepare_ligands()

    def prepare_receptors(self):
        """
        This function prepare receptors for plants 
        """
        receptors_prepared = []
        for rec in self.receptors:
            if not self._check_format_(rec):
                continue
            if self.verbose:
                print(f"\n\npreparing receptor {rec.name}")
            out = self.prepared_receptors_folder / f"{rec.stem}.mol2"
            self.run_spores(rec, out, mode="complete")
            receptors_prepared.append(out)
        self.receptors = receptors_prepared

    def prepare_ligands(self):
        """
        This function prepare ligands for plants 
        """
        ligands_prepared = []
        for lig in self.ligands:
            # .pdb .mol2 or .pdbqt
            if not self._check_format_(lig):
                continue
            if lig.suffix == ".mol2":
                self.fix_multiplefragment(lig)
            if self.verbose:
                print(f"\n\npreparing ligand {lig.name}")
            out = self.prepared_ligand_folder / f"{lig.stem}.mol2"
            self.run_spores(lig, out, mode="complete")
            ligands_prepared.append(out)
        self.ligands = ligands_prepared

    def get_docking_executable(self, docking_program):

        docking_engine = {
            "spores": "spores",  # linux & winc
            "plants": "plants",  # linux & win
        }

        exe_name = docking_engine[docking_program]
        return self._get_executable_(exe_name)

    def run_plants(self, out_path, lig_path, rec_path):

        conf_file = self.out_folder / rec_path.stem / \
            f"{rec_path.stem}-{lig_path.stem}-conf.txt"

        keywords = {
            "protein_file": rec_path,
            "ligand_file": lig_path,
            "output_dir": out_path,
        }

        # print(reader.read_conf(self.conf, "plants"))
        plants_keywords = reader.read_conf(self.conf, self.docking_program)
        plants_keywords = {**plants_keywords, **keywords}
        writer.write_conf(plants_keywords, "plants", conf_file)
        docking_executable = self.get_docking_executable(self.docking_program)
        command = f"{docking_executable} --mode screen {conf_file}"
        os.system(command)
        # os.popen(command).read()

    def run_screening(self, init_screening=True):

        if init_screening:
            self.init_screening()

        if self.prepare:
            self.logger("#Preparing ligand and receptors")
            self.prepare_screening()

        self.logger("#" * 72)
        self.logger(f"#Running virtual screening using {self.docking_program}")
        self.logger("#" * 72)

        for rec in self.receptors:
            if rec.suffix != ".mol2":
                continue
            list_args = []
            for lig in self.ligands:
                if lig.suffix != ".mol2":
                    continue
                if self.verbose:
                    self.logger(f"#Docking {lig.name} in {rec.name}")
                #                if os.listdir(self.ligand_folder).index(lig) <= x:
                #                    continue
                out_path = self.out_folder / rec.stem / lig.stem
                # os.mkdir(out_path)
                t0 = time()

                # self.run_plants(out_path, lig, rec)
                list_args.append((out_path, lig, rec))
                self.logger(
                    f"{self.docking_program},{rec.stem},{lig.stem},{time()-t0:.2f}")

            job_input_dock_lig = tuple(
                [tuple([self.run_plants, args]) for args in list_args]
            )
            multi_threading(job_input_dock_lig, 0, self.run_multiprocess)

    def run_multiprocess(self, func, args):
        return func(*args)

    def run_spores(self, input_file, output_file=None, mode="complete"):
        """
        Run a spores instance

        Parameters
        ----------
        input_file : pathlib.Path
            DESCRIPTION.
        output_file : pathlib.Path, optional
            DESCRIPTION. The default is None.
        mode : str, optional
            Spores mode.VAilbles [completemol2,
                                  reprot,
                                  settypes,
                                  readbonds,
                                  protstates,
                                  stereo,
                                  tautomers,
                                  ketoenol,
                                  ringconfs],
            The default is 'settypes'.

        Returns
        -------
        None.

        """
        if input_file == output_file:
            output_file = None

        if not output_file:
            output_file = input_file.parent / f"{input_file.stem}-prep.mol2"

        if sys.platform == "linux" or sys.platform == "linux2":
            mode = "complete"
        elif sys.platform == "win32":
            if input_file.suffix == ".pdb":
                mode = "completepdb"
            elif input_file.suffix == ".mol2":
                mode = "completemol2"

        docking_executable = self.get_docking_executable("spores")
        command = f"{docking_executable} --mode {mode} {input_file} {output_file}"

        # os.popen(command).read()
        return os.system(command)


class LedockScreening(Screening):
    def prepare_screening(self):
        """
        Look into receptors folder and ligands folder and prepare each
        file for docking. This use a SPORE for preparing ligand and receptor

        Parameters
        ----------
        verbose : bool, optional
            Print what its doing. The default is False.

        Returns
        -------
        None.

        """

        self.prepared_folder(folder_name="ledock")
        self.prepare_receptors()
        self.prepare_ligands()

    def prepare_receptors(self):
        """
        This function prepares receptors for ledock
        """
        receptors_prepared = []
        for rec in self.receptors:
            # pdb,mol2,pdbq,pdbqs,pdbqt, possibly pqr,cif
            receptors_prepared = []
            out = self.prepared_receptors_folder / f"{rec.stem}.pdb"
            self.run_lepro(rec)
            copy("./pro.pdb", out)
            os.remove("./pro.pdb")
            receptors_prepared.append(out)
        self.receptors = receptors_prepared

    def prepare_ligands(self):
        """
        This function prepares ligands for ledock
        """
        ligands_prepared = []
        for lig in self.ligands:
            # .pdb .mol2 or .pdbqt
            self.fix_multiplefragment(lig)
            out = self.prepared_ligand_folder / \
                f"{lig.stem}-prepared_ledock.mol2"
            self.run_spores(lig, out, mode="settypes")
            ligands_prepared.append(out)

        self.write_ligand_list()
        self.ligands = self.ligands_prepared

    def write_ligand_list(self):
        self.ligand_list = self.out_folder / "ligand.list"
        with open(self.ligand_list, "w") as f:
            for lig in self.ligands_prepared:
                f.write(f"{lig}\n")
        return self.ligand_list

    def get_docking_executable(self, docking_program):
        docking_engine = {
            "spores": "spores",  # linux & winc
            "lepro": "lepro",  # linux & win
            "ledock": "ledock",  # linux & win
        }

        exe_name = docking_engine[docking_program]
        return self._get_executable_(exe_name)

    def run_ledock(self, rec_path):

        rec_out = self.out_folder / rec_path.stem
        keywords = {
            "Receptor": str(rec_path.relative_to(rec_out)),
            "Ligands list": str(self.ligand_list.relative_to(rec_out)),
        }

        ledock_keywords = reader.read_conf(self.conf, "ledock")
        ledock_keywords = {**ledock_keywords, **keywords}

        conf_file = rec_out / f"{rec_path.stem}-conf.txt"
        writer.write_conf(ledock_keywords, "ledock", conf_file)

        ledock = self.get_docking_executable("ledock")
        os.chdir(rec_out)

        command = f"{ledock} {conf_file.name}"
        # command = f'{ledock} {conf_file.resolve()}'

        # os.popen(command).read()
        os.system(command)
        return True

    def run_screening(self, init_screening=True):

        if init_screening:
            self.init_screening()

        if self.prepare:
            self.logger("#Preparing ligand and receptors")
            self.prepare_screening()

        if self.verbose:
            self.logger("#" * 72)
            self.logger(
                f"#Running virtual screening using {self.docking_program}")
            self.logger("#" * 72)

        for rec in self.receptors:
            if rec.suffix != ".pdb":
                continue
            t0 = time()
            self.run_ledock(rec)
            self.logger(
                f"{self.docking_program},{rec.stem},null,{time()-t0:.2f}")

    def run_lepro(self, input_file, flag=""):
        """
        Run a spores instance

        Parameters
        ----------
        input_file : pathlib.Path
            DESCRIPTION.
        output_file : pathlib.Path, optional
            DESCRIPTION. The default is None.
        flag : str, optional
            -rot  [[chain] resid] align principal axes of 
                    the binding site with Cartesian                                           
            -metal keep ZN/MN/CA/MG
            -metal -p redistribute metal charge to protein 

        Returns
        -------
        None.

        """
        if input_file.suffix != ".pdb":
            print(f"{input_file.name} has not supported formats")
            print("lepro only support pdb format")
            # raise ValueError('lepro only support pdb format')
            return None

        lepro = self.get_docking_executable("lepro")
        command = f"{lepro} {input_file} {flag}"
        # os.popen(command).read()
        os.system(command)
        return 1

    def run_spores(self, input_file, output_file=None, mode="settypes"):
        """
        Run a spores instance

        Parameters
        ----------
        input_file : pathlib.Path
            DESCRIPTION.
        output_file : pathlib.Path, optional
            DESCRIPTION. The default is None..
        mode : str, optional
            Spores mode.VAilbles [completemol2,
                                  reprot,
                                  settypes,
                                  readbonds,
                                  protstates,
                                  stereo,
                                  tautomers,
                                  ketoenol,
                                  ringconfs],
            The default is 'settypes'.

        Returns
        -------
        None.

        """
        if not output_file:
            output_file = input_file.parent / f"{input_file.stem}-prep.mol2"
        docking_executable = self.get_docking_executable("spores")
        command = f"{docking_executable} --mode {mode} {input_file} {output_file}"
        os.popen(command).read()
        # os.system(command)
        return 1


vsprotocol = ScreeningFactory()
vsprotocol.register_program("vina", VinaScreening)
vsprotocol.register_program("fwavina", VinaScreening)
vsprotocol.register_program("gwovina", VinaScreening)
vsprotocol.register_program("psovina", VinaScreening)
vsprotocol.register_program("qvina2", VinaScreening)
vsprotocol.register_program("qvina-w", VinaScreening)
vsprotocol.register_program("smina", VinaScreening)

vsprotocol.register_program("plants", PlantsScreening)
vsprotocol.register_program("ledock", LedockScreening)
vsprotocol.register_program("dock", DOCKScreening)

# %%
if __name__ == "__main__":
    # =============================================================================
    # Test Vina
    # =============================================================================
    class Test_VinaScreen:

        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        out = Path("../../data/out")
        conf = Path("../../data/conf.txt")
        docking_program = "vina"

        def test_prepared_folder(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program)
            s.prepared_folder()

        def test_prepare_receptors(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program)
            s.prepared_folder()
            s.prepare_receptors()

        def prepare_ligands(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program)
            s.prepared_folder()
            s.prepare_ligands()

        def prepare_screening(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program)
            s.prepare_screening()

        def test_get_docking_executable(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program)
            docking_engine = {
                "fwavina",
                "vina",
                "gwovina",
                "ledock",
                "psovina",
                "qvina2",
                "qvina-w",
                "smina",
            }
            for docking_program in docking_engine:
                s.get_docking_executable(docking_program)

        def test_run_vina(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program)

            s.run_vina(self, "out_path", "lig_path", "rec_path")

        def test_run_screening_prepare_true(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program, prepare=True)
            s.run_screening()

        def test_run_screening_prepare_False(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.docking_program, prepare=False)
            s.run_screening()

    def test_vina():
        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        # prepared_ligands = Path("../../data/prepared/prepared_ligands_vina")
        # prepared_receptors = Path(
        # "../../data/prepared/prepared_receptors_vina")
        out = Path("../../data/out")
        conf = Path("../../data/config_vina_ex1.txt")
        docking_program = "vina"
        s = VinaScreening(
            ligands,
            receptors,
            out,
            conf,
            docking_program,
            prepare=True,
            # prepared_ligand_folder=prepared_ligands,
            # prepared_receptors_folder=prepared_receptors,
            exe_file="../exe.paths",
        )
        # s.init_screening()
        # s.prepared_folder('vina')
        # s.prepare_ligands()
        # s.prepare_receptors()
        # s.prepare_screening()
        # s.prepare_screening(verbose=1)
        s.get_docking_executable("vina")
        s.run_screening()

    # test_vina()
    #
    # %% =============================================================================
    # Plants test
    # =============================================================================
    def test_plants():
        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        # ligands = Path("../../data/prepare/prepared_ligands_plants")
        # receptors = Path("../../data/prepare/prepared_receptors_plants")
        out = Path("../../data/our-plants")
        conf = Path("../../data/config_plants_speed4.txt")
        docking_program = "plants"
        s = PlantsScreening(
            ligands, receptors, out, conf=conf, docking_program=docking_program, verbose=True, prepare=True, exe_file="../exe.paths"
        )
        # s.init_screening()
        # s.prepared_folder('plants')
        # s.prepare_ligands()
        # s.prepare_receptors()
        s.prepare_screening()
        # i = [l for l in ligands.iterdir()][0]
        # o = i.parent / f"{i.stem}-prep.mol2"
        # s.run_spores(i, o)
        # s.get_docking_executable("plants")
        s.run_screening(init_screening=True)

    # test_plants() 

    # %%===========================================================================
    # Ledock Test
    # =============================================================================
    # def test_ledock():
    #     ligands = Path("../../data/ligands")
    #     receptors = Path("../../data/receptor")
    #     prepared_ligands = Path("../../data/prepared/prepared_ligands_ledock")
    #     prepared_receptors = Path("../../data/prepared/prepared_receptors_ledock")
    #     out = Path("../../data/out")
    #     conf = Path("../../data/config_ledock_sample.txt")
    #     docking_program = "ledock"
    #     s = LedockScreening(
    #         ligands,
    #         receptors,
    #         out,
    #         conf,
    #         docking_program,
    #         # prepared_ligand_folder=prepared_ligands,
    #         # prepared_receptors_folder=prepared_receptors,
    #     )
    #     s.init_screening()
    #     s.prepare_screening(verbose=1)
    #     s.get_docking_executable("ledock")
    #     s.run_screening()

    # %%

    # ligands = Path('../../data/ligands')
    # receptors = Path('../../data/receptor')
    # prepared_ligand_folder = Path('../../data/prepared/prepared_ligands_vina')
    # prepared_receptors_folder = Path('../../data/prepared/prepared_receptors_vina')
    # out = Path('../../data/out-test-plants-multiprocess')
    # conf = Path('../../data/config_plants_speed4.txt')
    # docking_program = 'plants'

    # screening = vsprotocol.get_program(docking_program)
    # s = screening(ligands,
    #               receptors,
    #               out,
    #               conf,
    #               docking_program,
    #               prepared_ligand_folder,
    #               prepared_receptors_folder)

    # s.init_screening()
    # # s.prepare_screening(verbose=1)
    # s.run_screening()
    # #%%
    # ligands = Path('../../data/ligands')
    # receptors = Path('../../data/receptor')
    # prepared_ligand_folder = Path('../../data/prepared/prepared_ligands_plants')
    # prepared_receptors_folder = Path('../../data/prepared/prepared_receptors_plants')
    # out = Path('../../data/out-vina')
    # conf = Path('../../data/conf.txt')
    # docking_program = 'plants'

    # screening = vsprotocol.get_program(docking_program)
    # s = screening(ligands,
    #               receptors,
    #               out,
    #               conf,
    #               docking_program,
    #               prepared_ligand_folder,
    #               prepared_receptors_folder)

    # s.init_screening()
    # # s.prepare_screening(verbose=1)
    # s.run_screening()

    # =============================================================================
    # Test DOCK
    # =============================================================================

    def test_DOCK():
        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        # ligands = Path("../../data/prepare/prepared_ligands_dock")
        # receptors = Path("../../data/prepare/prepared_receptors_dock")
        out = Path("../../data/out2-dock")
        conf = Path("../../data/config_dock_flex_sample.txt")
        docking_program = "dock"
        s = DOCKScreening(
            ligands, receptors, out, conf=conf, docking_program=docking_program, verbose=True, prepare=True, exe_file="../exe.paths"
        )
        # s.init_screening()
        # s.prepared_folder('DOCK')
        # s.prepare_ligands()
        # s.prepare_receptors()
        # s.prepare_screening()
        # i = [l for l in ligands.iterdir()][0]
        # o = i.parent / f"{i.stem}-prep.mol2"
        # s.run_spores(i, o)
        # s.get_docking_executable("DOCK")
        s.run_screening()

    # test_DOCK() 