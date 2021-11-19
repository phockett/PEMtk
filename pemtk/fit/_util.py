
# Get and set generic class args and/or local fn values
# IN PROGRESS....
# ALSO, have done this before....?  ePSproc wavefn plotting case?
# ALSO, use a class factory here!!!

import numpy as np

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


# TODO: revise and wrap this properly!
# Should push corrected phase back to original table as new dim?
# Phase correction/shift function

def phaseCorrection(dfWide, dfRef = None, refParam = None, wrapFlag = True, phaseLabel = 'p', absFlag = False):  #, drop_level=False, renameLevel=True):
    """
    Phase correction/shift/wrap function.

    Prototype from test code:

    - Assumes full Pandas tabulated wide-form dataset as input.
    - Supply dfRef to use reference phase (abs phase values), otherwise will be relative with refParam set to zero.
    - wrapFlag: wrap to -pi:pi range? Default True.
    - absFlag: set abs values (drop signs)? Default False, otherwise sets all values to abs().

    TODO: implement dim preservation? Currently handled by calling fn., and returned values will have Type dim dropped here.
            Setting "drop_level=False" to xs() would fix this.
            Stated to implement, but skipped for now.

    """

    phasesIn = dfWide.xs(phaseLabel,level='Type')  # Set phase data
    phaseCorr = phasesIn.copy()

#     print(refParam)

    if refParam is None:
        refParam = phasesIn.columns[0]  # Default ref phase

    refPhase = dfWide[refParam].xs(phaseLabel,level='Type')  #,drop_level=drop_level)
#     print(refPhase)

    # For abs ref phase, set that too
#     absFlag = False
    if dfRef is not None:
        refPhase = refPhase - dfRef[dfRef['Param'] == refParam].xs(phaseLabel,level='Type')['value'].item()
#         absFlag = True

    print(f"Set ref param = {refParam}")

    # Subtract (shift) by refPhase
    phaseCorr = phaseCorr.subtract(refPhase, axis='index')  # Subtract ref phase, might be messing up sign here?

    # Rename dim level? Only valid if level is NOT dropped
    # if renameLevel and ~drop_level:
    #     phaseCorr.index = dfOut.index.set_levels(['m','pc'], level = 'Type')  # Set 'pc' Type - NOTE that index keeps all types.

    # Rectify phases...? Defined here for -pi:pi range.
    if wrapFlag:
        # phaseCorrRec = phaseCorr.apply(lambda x: np.sign(x)*np.mod(np.abs(x),np.pi)) # This will wrap towards zero, should be OK for zero ref phase case.

        # Use arctan, defined for -pi:pi range
        phaseCorr =  np.arctan2(np.sin(phaseCorr), np.cos(phaseCorr))

    # Set abs values? (Do this last!)
    if absFlag:
        phaseCorr = np.abs(phaseCorr)


    return phaseCorr
