# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/bin/runAdt.py,v 1.33 2008/07/30 19:02:47 vareille Exp $
# $Id: runAdt.py,v 1.33 2008/07/30 19:02:47 vareille Exp $

# AutoDockTools can be launched from a python shell like this:
#import AutoDockTools; AutoDockTools.runADT()

import sys
import AutoDockTools

if '__IP' in dir(): # ipython
    ownInterpreter = False
else:
    ownInterpreter = True 

if '__file__' in locals():
    AutoDockTools.runADT(sys.argv, ownInterpreter=ownInterpreter, AdtScriptPath=__file__)
else:
    AutoDockTools.runADT(sys.argv, ownInterpreter=ownInterpreter)
