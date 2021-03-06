#
# Author Michel F. Sanner  (may 2001)          Copyright M. Sanner, TSRI
#
# $Id: callback.py,v 1.10.12.2 2017/04/06 23:02:57 annao Exp $
#
# $Author: annao $
#
import traceback
import types

class CallbackManager:
    """Class to manage a list of callback functions"""

    def __init__(self):
	self.callbacks = []


    def FindFunctionByName(self, funcName):
	"""find a function with a given name in a list of functions"""

        for f in self.callbacks:
	    if f.__name__==funcName: return f
	return None


    def SetCallback(self, func):
	"""Delete all and set a callback fuction"""

	assert func is None or callable(func)
        if func is None:
            self.callbacks = []
        else:
            self.callbacks = [func, ]


    def AddCallback(self, func):
	"""Add a callback fuction"""

	assert callable(func)
	self.callbacks.append(func)


    def CallCallbacks(self, *args, **kw):
	"""call all callback fuctions"""

        results = []
	for func in self.callbacks:
            try:
                results.append( func(*args, **kw) )
                
            except:
                print('ERROR ************************************************')
                traceback.print_exc()
                return 'ERROR'

        return results

            
    def ListCallbacks(self):
        for func in self.callbacks:
            print(func.__name__,func)


    def RemoveCallback(self, func):
	"""Delete a callback fuction"""
        if type(func)==bytes:
            func = self.FindFunctionByName(func)
            if func is None: return "function %s not found"%func
        if func in self.callbacks:
            self.callbacks.remove(func)
        else:
            return "function %s not found"%func.__name__


class CallbackFunction:
    """Class to allow to specify arguments to a callback function"""

    def __init__(self, function, *args, **kw):
        self.function = function
        if hasattr(function, '__name__'):
            self.__name__ = function.__name__
        elif hasattr(function, 'name'):
            self.__name__ = function.name
        elif hasattr(function, '__class__') and hasattr(function.__class__, '__name__'):
            self.__name__ = f.__class__.__name__
        else:
            self.__name__ = 'noname' 
        self.args = args
        self.kw = kw

    def __call__(self, *args, **kw):
        args = self.args + args
        kw.update(self.kw)
        return self.function(*args, **kw)

CallBackFunction = CallbackFunction
