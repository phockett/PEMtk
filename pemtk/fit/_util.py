
# Get and set generic class args and/or local fn values
# IN PROGRESS....
# ALSO, have done this before....?  ePSproc wavefn plotting case?
# ALSO, use a class factory here!!!

# Set class args at init
def setClassArgs(self,args):
    for argName in args:
        if (argName!='self') and (argName!='kwargs'):
            setattr(self, argName, args[argName])
        elif argName=='kwargs':
            for k,v in args[argName].items():
                setattr(self, k, v)

# Check class args and overwrite if passed...?
# def checkClassArgs(self,args):
#     for argName in args:
#         if (argName!='self') and (argName!='kwargs'):
#             if args[argName] is not None:
#                 setattr(self, argName, args[argName])
#         elif argName=='kwargs':
#             for k,v in args[argName].items():
#                 if args[argName] is not None:
#                     setattr(self, k, v)



def lmmuListStrReformat(lmmuList):
    """Convert list of tuple labels to short str format"""
    # Parameter names [a-z_][a-z0-9_]*, so replace - sign with
    return ['_'.join(str(ele).replace('-','n') for ele in sub) for sub in lmmuList]


# Set default params for fitting analysis
# May want to change to use set/get style? Or use a decorator?
def _setDefaultFits(self, dataRange):

    # Set default indexes
    if dataRange is None:
        dataRange = [0, self.fitInd]

    return dataRange
