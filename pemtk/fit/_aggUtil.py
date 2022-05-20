# Routines for dealing with aggregate matrix elements
#
# 17/04/22  v1  Developed for manuscript plotting, see http://stimpy:8966/lab/tree/work/pemtk_fitting_runs_April2022/analysis_dev/MFrecon_manuscript_fig_generation_170422-Stimpy-DEV.ipynb

import numpy as np

# Functional form with types
from epsproc.util.listFuncs import dataTypesList

def setAggCompData(data, keys = {'comp':['m','p'],'compC':['n','pc']}):
    """
    Set complex forms for aggregate results from [mag,phase] columns.

    Parameters
    ----------
    data : pd.dataframe
        Data to convert.

    keys : dict, optional, default = {'comp':['m','p'],'compC':['n','pc']}
        Dict of keys for output, and [mag,phase] columns to convert.

    TODO: generalise further?

    """

    # data['comp'] = data['m']*np.exp(data['p']*1j)
    # data['compC'] = data['m']*np.exp(data['pc']*1j)

    for k,v in keys:
        data[k] = data[v[0]]*np.exp(data[v[1]]*1j)

    return data


def setAggMatE(self, simpleForm = False):
    """Set aggregate results to matE format"""

    # Final results & reformat
    meanMatE = self.paramsSummary['agg'].query("(Agg == 'mean')")

    # Replace index...
    meanMatE.index = self.lmmu['Index']

    meanMatE = setAggCompData(meanMatE)  # Set complex forms

    if simpleForm:
        # Relabel, paper format - see ._util.lmmuListStrReformat()
        labels = meanMatE.droplevel(['Cont','Targ','Total','mu','it']).index
        lmmuList = labels
        strLabels = [','.join(str(ele) for ele in sub) for sub in lmmuList]

        meanMatE = meanMatE.assign(labels=strLabels)

        meanMatE = meanMatE.droplevel(['Targ','Total','mu','it'])  #.reset_index()  #.set_index('labels')

    # Set output
    self.data['agg']['matEpd'] = meanMatE


def aggToXR(self, key = 'agg', cols=['comp','compC'], EkeList = [1.1], dType = 'matE',
               refKey = None, returnType = 'ds'):
    """
    Pull columns from PD dataframe & stack to XR dataset.


    TODO:
    - EkeList from input data subset?
    - Use existing routines for more flexible dim handling? E.g. pemtk.sym._util.toePSproc
    - More returnTypes, currently set for single dataset or set of arrays (per col)

    """

    refDims = dataTypesList()[dType]['def']

    # Get mean matE if not set
    if not f'{dType}pd' in self.data[key].keys():
        self.setAggMatE()
        # setAggMatE(data)
#         setAggCompData


    dataSet = {}
    for col in cols:
#         print(col)
        # Basic XR array
        xrRaw = xr.DataArray(self.data[key][f'{dType}pd'][[col]]).unstack()
#         print(xrRaw)

        # Reformat
        dataSet[col] = xrRaw.stack(refDims(sType='sDict')).squeeze(drop=True).expand_dims({'Eke':EkeList})

    if refKey is not None:
        dataSet[refKey] = self.data[refKey][dType]

    # Stack to DataSet
#     return xr.Dataset(dataSet)
    if returnType == 'ds':
        self.data[key][dType] = xr.Dataset(dataSet)

        if self.verbose:
                print(f"Set XR dataset for self.data['{key}']['{dType}']")

    elif returnType == 'da':
        for k,v in dataSet.items():
            label = f'{key}_{k}'
            self.data[label] = {dType:v}

            if self.verbose:
                print(f"Set XR dataarray for self.data['{label}']['{dType}']")
