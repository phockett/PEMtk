# Routines for dealing with aggregate matrix elements
#
# 17/07/22  v2  Debugged plus a few extensions.
# 17/04/22  v1  Developed for manuscript plotting, see http://stimpy:8966/lab/tree/work/pemtk_fitting_runs_April2022/analysis_dev/MFrecon_manuscript_fig_generation_170422-Stimpy-DEV.ipynb

import numpy as np
import xarray as xr

# Functional form with types
from epsproc.util.listFuncs import dataTypesList

# Safe restacker
from epsproc.util.misc import restack


def setAggCompData(dataIn, keys = {'comp':['m','p'],'compC':['n','pc']}):
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

    data = dataIn.copy()   # Use .copy() to fix slice warning in some PD versions (this will always be a new col)
                           # Using .assign() might be better?

    for k,v in keys.items():
        data[k] = data[v[0]]*np.exp(data[v[1]]*1j)

    return data


def setAggMatE(self, key = 'agg', dataOut = None,
                compDataLabels = {'comp':['m','p'],'compC':['n','pc']},
                simpleForm = False, dropLabelsList = ['Cont','Targ','Total','mu','it'],
                dropLevelsList = ['Targ','Total','it']):
    """
    Set aggregate results to matE format (Pandas)

    If key='ref' use self.data[self.subKey]['matE'] instead of aggregate data

    18/07/22 - quickly hacked in ref data case for consistent results tabulations, probably already have this stuff elsewhere.
                See also pemtk.fit._conv.pdConvRef() and self.setMatEFit()
    """
    if dataOut is None:
        dataOut = 'matEpd'

    if key is not 'ref':
        # Final results from aggregate data & reformat


        meanMatE = self.paramsSummary[key].query("(Agg == 'mean')")

        # Replace index...
        meanMatE.index = self.lmmu['Index']

        meanMatE = setAggCompData(meanMatE, compDataLabels)  # Set complex forms

    else:
        # Set as ref. data & reformat
        # 18/07/22 - quickly hacked this in for consistent results tabulations.
        # TODO: probably already have this elsewhere? Should work from self.dfRef?
        # if dataOut is None:
        #     dataOut = 'matEpdRef'

        meanMatE = self.pdConvSetFit(matE = self.data[self.subKey]['matE'])

        # Set cols for complex & abs/phase forms.
        meanMatE.columns = ['comp']
        meanMatE[['m']] = np.abs(meanMatE[['comp']])
        meanMatE['p'] = np.angle(meanMatE[['comp']])


    if simpleForm:
        # Relabel, paper format - see ._util.lmmuListStrReformat()
        # labels = meanMatE.droplevel(['Cont','Targ','Total','mu','it']).index
        # labels = meanMatE.droplevel(['Cont','Targ','Total','it']).index   # keep mu
        labels = meanMatE.droplevel(dropLabelsList).index   # With passed args
        lmmuList = labels
        strLabels = [','.join(str(ele) for ele in sub) for sub in lmmuList]

        meanMatE = meanMatE.assign(labels=strLabels)

        # meanMatE = meanMatE.droplevel(['Targ','Total','mu','it'])  #.reset_index()  #.set_index('labels')
        # meanMatE = meanMatE.droplevel(['Targ','Total','it'])  # Keep mu
        meanMatE = meanMatE.droplevel(dropLevelsList)  # With passed args.

    # Set output
    if not key in self.data.keys():
        self.data[key] = {}

    self.data[key][dataOut] = meanMatE

    if self.verbose['main']:
        if key is not 'ref':
            print(f'Set reformatted aggregate data to self.data[{key}][{dataOut}].')
        else:
            print(f'Set reformatted ref data to self.data[{key}][{dataOut}].')


def aggToXR(self, key = 'agg', cols = {'comp':['m','p'],'compC':['n','pc']},   # cols=['comp','compC'],
            EkeList = [1.1], dType = 'matE', conformDims=True,
            refKey = None, returnType = 'ds', simpleForm = False):
    """
    Pull columns from PD dataframe & stack to XR dataset.


    cols : dict, optional, default = {'comp':['m','p'],'compC':['n','pc']}
        Dict of keys for output items/columns, and [mag,phase] columns to convert.


    TODO:
    - EkeList from input data subset?
    - Use existing routines for more flexible dim handling? E.g. pemtk.sym._util.toePSproc
    - More returnTypes, currently set for single dataset or set of arrays (per col)

    UPDATE 19/07/22 - implemented ep.misc.restack(), which includes dim checking and expansions.
                    Set `conformDims` True/False for the latter.
                    Eke dim still handled separately here.
                    NOTE: conformDims=False with refKey only works reliably for da return type, otherwise may fail at dataset stacking.

    """

    # refDims = dataTypesList()[dType]['def']

    # Get mean matE if not set
    if not f'{dType}pd' in self.data[key].keys():
        self.setAggMatE(keys = cols, simpleForm = simpleForm)
        # setAggMatE(data)
#         setAggCompData


    dataSet = {}
    dimSets = {}
    for col in cols.keys():
#         print(col)
        # Basic XR array
        xrRaw = xr.DataArray(self.data[key][f'{dType}pd'][[col]]).unstack()
#         print(xrRaw)

        # Reformat
        # dataSet[col] = xrRaw.stack(refDims(sType='sDict')).squeeze(drop=True).expand_dims({'Eke':EkeList})

        # UPDATE 19/07/22 - use ep.misc.restack(), which includes dim checking and expansions.
        # TODO: pass/use dim maps as returned here too?
        #   dataSet[col], dims = restack(xrRaw, refDims = refDims)  # OK, with missing dims
        #   dataSet[col], dims = restack(xrRaw, refDims = refDims, conformDims=True)  # OK, but not adding dims?
        #   dataSet[col], dims = restack(xrRaw, refDims = 'matE', conformDims=True)  # OK, adds dims
        dataSet[col], dimSets[col] = restack(xrRaw, refDims = dType, conformDims=conformDims)  # OK, adds dims with flag.

    if refKey is not None:
        dataSet[refKey] = self.data[refKey][dType]

    # Stack to DataSet
#     return xr.Dataset(dataSet)

    # Set output
    if not key in self.data.keys():
        self.data[key] = {}

    if returnType == 'ds':
        try:
            self.data[key][dType] = xr.Dataset(dataSet)

        except Exception as e:
            print(f"*** Caught exception at xr.Dataset stacking: {e}")
            print("Try returnType = 'da' to avoid this and handle dims manually.")
            print("Returning dataSet dictionary instead...")

            return dataSet

        if self.verbose:
                print(f"Set XR dataset for self.data['{key}']['{dType}']")

    elif returnType == 'da':
        for k,v in dataSet.items():
            label = f'{key}_{k}'
            self.data[label] = {dType:v}

            if self.verbose:
                print(f"Set XR dataarray for self.data['{label}']['{dType}']")
