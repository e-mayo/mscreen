# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 17:35:34 2020

Prepare receptor from python 3 using the python 2.7 interpreter of ADFRsuite

@author: emayo
"""
from pathlib import Path
import os

ADFRsuite_path = "C:\Program Files (x86)\ADFRsuite-1.0"

ADFRsuite_path = Path(ADFRsuite_path)
python27 = ADFRsuite_path / "python.exe"
utils = ADFRsuite_path / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24"
prepare_receptor4 = utils / "prepare_receptor4.py"
#%%
receptor = Path("receptor_pdb3b2x.ent_out_ok.pdb")

print(os.popen(f""""{str(python27)}" "{prepare_receptor4}" -r {receptor} -o {receptor.stem}.pdbqt -v -p ZN""").read())

