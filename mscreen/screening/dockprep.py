import chimera
from DockPrep import prep
import sys
outfilename = sys.argv[4]

models = chimera.openModels.list(modelTypes=[chimera.Molecule])
prep(models,
    addCharges=True)


# if outfilename.split('.')[-1] == 'mol2'
from WriteMol2 import writeMol2
writeMol2(models, outfilename)



