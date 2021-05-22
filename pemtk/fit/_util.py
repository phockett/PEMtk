
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
