# PEMtk, converters for fit results > Pandas
#
# 15/06/21  v1


import pandas as pd
import numpy as np
from epsproc import multiDimXrToPD

# TODO: should be able to simplify & automate with results.__dict__  - May not be comprehensive?
# NOW REINDEX BY FIT # & Type, this makes sense for wide <> long conversions
# For stacked case, set as tuple dictionary index
def pdConv(self, fitVars = ['success', 'chisqr', 'redchi'], paramVars = ['value', 'stderr', 'vary', 'expr'],
            dataRange = None, batches = None):
    """
    Basic conversion for set of fit results > Pandas, long format.

    Extract fit and parameter results from lmFit objects and stack to PD dataframe.


    Parameters
    ----------

    fitVars : optional, list, default = ['success', 'chisqr', 'redchi']
        Values to extract from lmfit result object (per fit).

    paramVars : optional, list, default = ['value', 'stderr', 'vary', 'expr']
        Values to extract from lmfit params object (per parameter per fit).

    dataRange : optional, list, default = None
        Range of indexes to use, defaults to [0, self.fitInd].

    batches : optional, int, default = None
        Additional batch of labelling for fits.
        - If int, label as ceil(fit #)/batches. E.g. batches = 100 will label fits per 100.
        - If list, use as labels per fit. (NOT YET IMPLEMENTED)

    Todo

    - Additional batching options, inc. by file for multiple read case.

    """

    # Set default indexes
    # if dataRange is None:
    #     dataRange = [0, self.fitInd]
    dataRange = self._setDefaultFits(dataRange)

    # Set vars
    dataDict = {}
    # variables = ['value', 'stderr']  #  Per parameter values.
    # fitVars = ['success', 'chisqr', 'redchi'] # These are per fit, not per param - should handle separately?

    outputIndex = ['Fit','Type','pn']  # Cols in output PD dataframe

    # Extract relevant data from lmfit params class objects & reindex
    for fitInd in range(dataRange[0], dataRange[1]):
        try:
            # Get fit vars
            fitDict = {k:getattr(self.data[fitInd]['results'], k) for k in fitVars}

            for n,i in enumerate(self.data[fitInd]['results'].params.items()):
        #     for n,i in enumerate(data.result.params.items()):
            #     print(n,i)

                pmType = i[0][0]
                dataDict[(fitInd, pmType, n)] = {j:getattr(i[1],j) for j in paramVars}
                dataDict[(fitInd, pmType, n)]['Param'] = i[0][2:]  # Use name + type for easier plotting later?
            #     dataDict[(fitInd, n)]['Type'] = i[0][0]
            #     dataDict[n]['Fit'] = fitInd  # As column

                dataDict[(fitInd, pmType, n)].update(fitDict)  # Add per fit items

        except KeyError as e:
            if self.verbose['sub']:
                print(f"*** Missing fit key {n}")

    # Stack to long-format PD
    dfLong = pd.DataFrame(dataDict).T
    # dfLong = pd.DataFrame.from_dict(dataDict).T  # Same result

    # Set index names
    dfLong.index.names = outputIndex

    # Set dType attrib for later.
    dfLong.attrs['dType'] = 'Params Long'

    # Set batches if specified
    if batches:
        if isinstance(batches,int):
            dfLong['batch'] = np.ceil((dfLong.index.get_level_values('Fit')+1)/(batches)).astype(int)  # Quick category for batch of fits
            # outputIndex.append('batch')

        # TODO: implement some other options here!

        else:
            print(f'** Batch type {type(batches)} not yet supported, skipping batching in pdConv() routine.')

    # Set ref values too, if present
    if hasattr(self,'params'):
        # for n,i in enumerate(self.params.items()):
        # #     print(n,i)
        #
        #     pmType = i[0][0]
        #     refDataDict[('ref', pmType, n)] = {j:getattr(i[1],j) for j in paramVars}
        #     refDataDict[('ref', pmType, n)]['Param'] = i[0][2:]  # Use name + type for easier plotting later?
        # #     dataDict[(fitInd, n)]['Type'] = i[0][0]
        # #     dataDict[n]['Fit'] = fitInd  # As column
        #
        # # Stack to long-format PD
        # dfRef = pd.DataFrame(refDataDict).T
        # dfRef.index.names = outputIndex
        # dfRef.attrs['dType'] = 'Params Ref'

        dfRef = self.pdConvRef(paramVars, outputIndex)  # Functionalised version

    else:
        dfRef = None
        print("Pandas reference table not set, missing self.params data.")


    return(dfLong, dfRef)


def pdConvRef(self, paramVars = ['value'], outputIndex = ['Fit','Type','pn']):
    """
    Convert reference params set to reference PD table.

    Basic routine stripped from main pdConv() method for reuse elsewhere.

    TODO: add flexibility here.

    """

    # If params are missing, try to set them
    if not hasattr(self,'params'):
        # Check if set in subset (only in code as of 29/11/21)
        try:
            self.params = self.data[self.subKey]['params']
            self.lmmu = self.data[self.subKey]['lmmu']

        except:
            print(f"self.params not set, setting ref values from self.data[{self.subKey}]['matE'] without constraints.")
            self.setMatEFit(paramsCons = {})

    refDataDict = {}

    for n,i in enumerate(self.params.items()):
    #     print(n,i)

        pmType = i[0][0]
        refDataDict[('ref', pmType, n)] = {j:getattr(i[1],j) for j in paramVars}
        refDataDict[('ref', pmType, n)]['Param'] = i[0][2:]  # Use name + type for easier plotting later?
    #     dataDict[(fitInd, n)]['Type'] = i[0][0]
    #     dataDict[n]['Fit'] = fitInd  # As column

    # Stack to long-format PD
    dfRef = pd.DataFrame(refDataDict).T
    dfRef.index.names = outputIndex
    dfRef.attrs['dType'] = 'Params Ref'

    return dfRef


def pdConvSetFit(self, matE, colDim = 'it'):
    """
    Restack matE to pd.DataFrame and force to 1D.

    Utility function for setting up fit parameter sets.

    """

    # Using PD conversion routine works, although may have issues with singleton dims again - should set suitable dummy dim here?
    # pdTest, _ = ep.multiDimXrToPD(testMatE, colDims='Eke', dropna=True)
    pdTest, _ = multiDimXrToPD(matE, colDims=colDim, dropna=True, squeeze = False)
    # pdTest, _ = ep.multiDimXrToPD(testMatE, colDims='Sym', dropna=True, squeeze = False)

    pdTest = pd.DataFrame(pdTest.stack(colDim))  # Stack to 1D format and force to DF

    return pdTest
