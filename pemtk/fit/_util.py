
# Get and set generic class args and/or local fn values
# IN PROGRESS....
# ALSO, have done this before....?  ePSproc wavefn plotting case?
# ALSO, use a class factory here!!!

import pandas as pd
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


# NOTE - this will only work to add a single extra level, and then will nest.
# More general solutions: https://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex
def addColLevel(df, newCol = 'ref', names = ['Dataset','Type']):
    """Add top-level column to Pandas dataframe (crude)"""

    oldCols = df.columns.to_list()
    newCols = pd.MultiIndex.from_product([[newCol], oldCols], names = names)
    df.columns = newCols


# Quick go at general remapper from dict - need to change by type (ind, col, val)? What about multiindex?
# Maybe with getattr(data.data['plots']['corrData'],'columns'), then by type? (Index or not)
# Should be a solved problem...?
def renameParams(data, mapDict, mapNames = ['lm'], mapType = 'col'):
    """
    Very basic column name reampper for Pandas DataFrame. Based on routine in SymHarm class.

    21/04/22    v1 for testing only.

    TODO: generalise this & consolidate!

    See also value remapping in paramPlot() routine, runs data.replace({'Param':self.lmmu[remap]}, inplace=True)

    """


    # With variable len remap
    newNames = []
    for k in data.columns:
#         print(k)

# If in attrs
#         if k in inputData.attrs['indexes'][names].keys():
#             newNames.append(inputData.attrs['indexes'][names][k])
#         else:
#             newNames.append(k)

# If passed
        if k in mapDict.keys():
            newNames.append(mapDict[k])
        else:
            newNames.append(k)

#     print(newNames)
#         test = self.coeffDF.copy()   # Display as is

    test = data.copy()   # With unstack and fillna
#     print(test.columns.rename(Params = newNames))
    test.columns = pd.MultiIndex.from_arrays([newNames], names = mapNames).get_level_values(0)  # Force to single level?
    # print(test)

    return test


# Basic renorm function, cf. _util.phaseCorrection()
def renormMagnitudes(dfWide):
    """
    Basic renormalisation of magnitudes so sum(mags**2) = 1

    Prototype from test code:

    - Assumes full Pandas tabulated wide-form dataset as input.
    - Renormed values appended to input dataframe, as Type=n

    TODO: implement dim preservation? Currently handled by calling fn., and assumes Type=m is present in index.

    """
    # Renorm magnitudes?
    attrsIn = dfWide.attrs.copy()

    # Note this currently only runs from dfWide data
    mSq = (dfWide**2).sum(axis=1)  # Set norm const as sum of squares (all Types)
    renormTest = (dfWide.xs('m',level='Type',drop_level=False)).divide(np.sqrt(mSq.xs('m',level='Type',drop_level=False)), axis=0)  # Divide out
    # Restack frames
    dfRenorm = pd.concat([dfWide, renormTest.rename({'m':'n'})]).sort_values(by=['Fit','Type'])  # Sort by Fits? Otherwise new vals appended to end

    # Propagate attrs
    dfRenorm.attrs = attrsIn
    dfRenorm.attrs['renorm'] = True

    return dfRenorm



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
    # For more control consider np.unwrap() and np.mod()
    # https://numpy.org/doc/stable/reference/generated/numpy.unwrap.html
    if wrapFlag:
        # phaseCorrRec = phaseCorr.apply(lambda x: np.sign(x)*np.mod(np.abs(x),np.pi)) # This will wrap towards zero, should be OK for zero ref phase case.

        # Use arctan, defined for -pi:pi range
        phaseCorr =  np.arctan2(np.sin(phaseCorr), np.cos(phaseCorr))

    # Set abs values? (Do this last!)
    if absFlag:
        phaseCorr = np.abs(phaseCorr)


    return phaseCorr
