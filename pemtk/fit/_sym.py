# PEMtk, routines for checking matrix element relations/symmetries and setting constraints
#
# 31/08/21, packaging from dev notebook v3 code. See ePS_matE_symmetrization_dev_280621.ipynb

#************* v3-pk 31/08/21, packaging v3 code as set of functions.

import pandas as pd
import numpy as np
from ._util import lmmuListStrReformat


# TODO: try simplifying sym definitions with functions. See https://lmfit.github.io/lmfit-py/constraints.html#advanced-usage-of-expressions-in-lmfit
# Quick test here, from ePS_matE_symmetrization_dev_280621.ipynb
# # Test adding function to lmfit params - note this requires params object to already be defined!
# def phaseWrap(x):
#     """Phase shift wrapper for -pi:pi"""
#
# #     return arctan2(sin(x), cos(x))          # This version not working in asteval run, although does work if set directly in params.
#     return np.arctan2(np.sin(x), np.cos(x))   # This version does run - namespace (local/global/asteval) issue?
#
# data.params._asteval.symtable['phaseWrap'] = phaseWrap


def symCheckDefns():
    """
    Set dictionary of lambda functions for matrix element symmetry relation checks.

    Currently set to check for identity, abs, phase & complex rotations (+/-pi/2).

    """

    # Nested version with additional fn string/constraint defn.
    # Here:
    #   - set 'transform' true/false to specify whether to operate on source data too.
    #   - Set constraint with 'm_' and 'p_' to match current params settings. Should unify this somehow!
    #   - NOTE: current version DOESN'T INCLUDE PHASE WRAPPING, so might cause issues with -pi:pi limits? Should test here with np.arctan2 solution as above.
    lams = { 'i':  {'name': 'identity',    # Restructure with nesting...? Might be overkill, but would be more flexible.
                    'lam': lambda x: x,
                    'transform': False,
                    'constraint': 'x'},
            'abs': {'name':'abs',
                    'lam': lambda x: np.abs(x),
                    'transform': True,
                    'constraint': 'm_x'},
            'phase': {'name':'phase',
                        'lam': lambda x: np.angle(x),
                        'transform': True,
                        'constraint': 'p_x'},
            'crot_p': {'name':'Complex rotation +pi/2',
                        'lam': lambda x: -x.imag + x.real*1j,
                        'transform': False,
    #                     'constraint': 'p_x + pi/2'},
                        'constraint': 'arctan2(sin(p_x+pi/2), cos(p_x+pi/2))'},   # Phase-wrapped version, seems to work OK in testing
            'crot_m': {'name':'Complex rotation -pi/2',
                        'lam': lambda x: x.imag - x.real*1j,
                        'transform': False,
    #                     'constraint': 'p_x - pi/2'},
                        'constraint': 'arctan2(sin(p_x-pi/2), cos(p_x-pi/2))'},   # Phase-wrapped version, seems to work OK in testing
           }

    return lams



def symCheck(self, pdTest = None, matE = None, colDim = 'it', lams = None, verbose = 1):
    """
    Check symmetrization of input matrix elements.

    Parameters
    ----------
    pdTest : pandas DataFrame, optional, default = None
        Matrix elements to check, as set in a Pandas table format.
        Currently expects 1D array of matrix elements, as set in setMatEFit().
        If None, matE will be used to create the test data.

    matE : Xarray, optional, default = None
        Matrix elements to check.
        If None, then uses self.data[self.subKey]['matE']

    colDims : dict, default = 'it'
        Quick hack to allow for restacking via ep.multiDimXrToPD, this will set to cols = 'it', then restack to 1D dataframe.
        This should always work for setting matE > fit parameters, but can be overridden if required.

    lams : dict, optional, default = None
        Dictionary of test lambda functions.
        If not set, will use lams = symCheckDefns()

    


    Returns
    -------
    dict
        Set of parameter mappings/constraints, suitable to use for self.setMatEFit(paramsCons = newDict)

    dict of DataFrames
        - 'unique', Reduced set of unique matrix elements only.
        - 'constraints', List of constraints (as per parameters dict).
        - 'tests', Full list of tests & relations found.


    TODO:

    - Wrap for class.
    - Input checks and set default cases (see setMatEFit()). Should tidy to single input and then check type?

    """

    # Quick and ugly wrap args for class - should tidy up here!
    if pdTest is None:

        if matE is None:
            # matE = self.subset
            matE = self.data[self.subKey]['matE']

        # pdTest, _ = multiDimXrToPD(matE, colDims=colDim, dropna=True, squeeze = False)
        # # pdTest, _ = ep.multiDimXrToPD(testMatE, colDims='Sym', dropna=True, squeeze = False)
        # pdTest = pd.DataFrame(pdTest.stack(colDim))  # Stack to 1D format and force to DF

        pdTest = self.pdConvSetFit(matE, colDim)  # Functional version.

    # Get lambda functions if not set
    if lams is None:
        lams = symCheckDefns()


    # Create empty frame to hold expressions
    dfExpr = pd.DataFrame().reindex_like(pdTest)

    for k,v in lams.items():
        dfExpr[k] = dfExpr[0].copy()

    # Try looping
    # Here:
    #   - Groupby [Cont,l] & loop over test items & conditions.
    #   - May be better/neater/faster ways with native pandas (corr, vector/matrix mult?), itertools, or numpy methods?
    #   - No way to pop a row item in pd?

    for df in pdTest.groupby(['Cont','l']):  # Treat (continua, l) independently
    #     print(df[0])  # Labels
    #     print(df[1])  # DF group
    #     print(df[1][0])  # DF group column
    #     print(df[1][0][1])  # Individual values/items

        for testItem in df[1][0].items():   # Items + labels
    #         print(testItem)
    #         print(pdTest1D.loc[testItem[0]])  # Can also use label as selector

    #         lmmuListStr = lmmuListStrReformat(df[1][0].index)  # Set full index
            testLabel = lmmuListStrReformat([testItem[0]])[0]  # Test item only
            df[1][0].drop(testItem[0], inplace=True)  # Pop test item to avoid self matches

            for k,v in lams.items():
    #             testVal = testItem[1].apply(v)
                testVal = v['lam'](testItem[1])

                if v['transform']:
                    dfTest = df[1][0].apply(v['lam']) == testVal
                else:
                    dfTest = df[1][0] == testVal

    #             # In some cases also want to change source df for testing
    #             if k in ['abs', 'phase']:
    #                 dfTest = df[1][0].apply(v) == testVal
    # #                 testExpr =
    # #             elif k.startswith('rot'):
    # #                 dfTest = df[1][0].apply(np.angle).round(phasePrecision) == testVal  # For rotation case need to compare phase val with rotated phase val.
    # #                 dfTest = df[1][0].apply(np.angle) == testVal
    #             else:
    #                 dfTest = df[1][0] == testVal

    #             dfTest.drop(test)

    #             print(dfTest)
    #             print(lmmuListStr[dfTest.values])
    #             if dfTest.sum() >0:  # Add sum test here, avoids issues with empty or missing values. ACTUALLY NOT REQUIRED
    #                 print(lmmuListStr[0])

                # Inplace setting
    #             dfExpr.where(~dfTest, testLabel, inplace = True)  # Not ~dfTest here, replaces FALSE values.

                # Version with col per test expr
    #             dfExpr[k].where(~dfTest, testLabel, inplace = True)

                # Try skipping reset for exprs already defined (e.g. abs)
    #             dfExpr[k][dfExpr[k].isna()].where(~dfTest, testLabel, inplace = True)  # Runs, but returns all NaN
    #             ind = dfExpr[k].isna()
                ind = dfExpr[k].isna()
    #             dfExpr[k].where(~(dfTest & ind), testLabel, inplace = True)  # OK, test case with basic relation label

                # With inspect + string replace for full fn sting output
                # This works nicely, but fn strings not quite as required for (abs,phase) fitting case!
    #             fnStr = str(inspect.getsourcelines(lams[k])[0][0]).split(':')[-1].strip().replace('x', testLabel).strip(',')

                fnStr = v['constraint'].replace('x',testLabel)
                dfExpr[k].where(~(dfTest & ind), fnStr, inplace = True)


    #**** Set final tables

    # dfCons = dfExpr['abs','phase'].copy()
    dfCons = pd.DataFrame(dfExpr['abs'])
    dfCons['phase'] = dfExpr['phase'].combine_first(dfExpr['crot_p']).combine_first(dfExpr['crot_m'])  # Combine phase terms with priority

    dfSymUn = pdTest[dfExpr['abs'].isna()]


    #**** Remap tables to dict to use for self.setMatEFit(paramsCons = newDict)
    pCons = dfCons.dropna().to_dict()  # This outputs 'abs' and 'phase' constraints dict in current form, just need to rename params in output and all is good.

    # pCons['abs']

    # Reformat with param labels (maybe easier to just set as a column above?)
    magDict = {f'm_{lmmuListStrReformat([k])[0]}': v for k,v in pCons['abs'].items()}
    phaseDict = {f'p_{lmmuListStrReformat([k])[0]}': v for k,v in pCons['phase'].items()}

    magDict.update(phaseDict)  # Combine


    #**** Returns

    return magDict, {'unique':dfSymUn, 'constraints':dfCons, 'tests':dfExpr}
