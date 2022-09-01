"""
PEMtk fitting base classes

Development for fitting class, take pemtk.dataClass as base and add some subselection & fitting methods.

10/05/21    v1b added fitting routines using ePSproc and lmfit libraries.
                Adapted from dev code notebook (functional forms), still needs some tidying up & wrapping for class.
                For dev code see http://127.0.0.1:8888/lab/tree/dev/PEMtk/fitting/fitting_routines_dev_PEMtk_300421.ipynb

13/10/20    v1

TODO:

- Clean up/finalise data scheme. Currently mix of dictionary style, self.data[dataset][datatype] and class attribs self.datatype. Better to go with former for flexibility, or latter for using as simple base class to wrap later? Former works with some existing plotting fns., but complicated.
- afblmMatEfit() defaults similarly messy, although basically working.
- More analysis tools for fitting & results to add, currently just self.fit() to run.

"""

import copy
from datetime import datetime as dt # Import datetime.datetime for now() function

import pandas as pd
import numpy as np
import xarray as xr

try:
    from lmfit import Parameters, Minimizer, report_fit

except ImportError:
    print("*** lmfit not found: data functions available, but not fitting.")


from pemtk.data.dataClasses import dataClass
from pemtk.util.env import isnotebook
from epsproc import matEleSelector, multiDimXrToPD, setADMs, setPolGeoms
from epsproc.geomFunc import afblmXprod, mfblmXprod
from epsproc.sphPlot import plotTypeSelector

# Set HTML output style for Xarray in notebooks (optional), may also depend on version of Jupyter notebook or lab, or Xr
# See http://xarray.pydata.org/en/stable/generated/xarray.set_options.html
# if isnotebook():
#     xr.set_options(display_style = 'html')

# def lmmuListStrReformat(lmmuList):
#     """Convert list of tuple labels to short str format"""
#     # Parameter names [a-z_][a-z0-9_]*, so replace - sign with
#     return ['_'.join(str(ele).replace('-','n') for ele in sub) for sub in lmmuList]

# Now moved to ._util, should move to ..util?
from ._util import lmmuListStrReformat
from ._stats import setPoissWeights

# Plotters - testing import routines here, include modules & flags
# from epsproc.plot import hvPlotters  # Should probably use this!
# from pemtk.util.hvPlotters import setPlotDefaults  # Pkg version, opts function only (as per tmo-dev)
# print(hvFlag)  # Fails
from pemtk.util import hvPlotters  # Imports full module, include setup routines.
# print(hvPlotters.hvFlag)  # OK
# print(hvFlag)  # Fails
# print(hvPlotters.hv)

try:
    import seaborn as sns
    snsFlag = True
except ImportError:
    print("*** Seaborn not found: some plotting functions may be unavailable.")
    snsFlag = False


class pemtkFit(dataClass):
    """
    Class prototype for pemtkFit class. Dev version builds on dataClass, and adds some basic subselection & fitting capabilities.
    """

    from ._aggUtil import aggToXR, setAggMatE
    from ._analysis import (analyseFits, fitHist, fitsReport, classifyFits, corrPlot, paramPlot, paramsReport, paramsCompare, paramFidelity,
                            phaseCorrection, _mergePFLong, _setData, _setWide)  #, scopeTest
    from ._conv import pdConv, pdConvRef, pdConvSetFit
    from ._filters import thresFits, _subsetFromXS, getFitInds
    from ._io import processedToHDF5, writeFitData, loadFitData, _pickleData, _writePDData, _writeXRData, setTimeStampedFileName
    from ._parallel import multiFit
    from ._plotters import BLMfitPlot, lmPlotFit, BLMsetPlot, hvSave
    # from ._stats import setPoissWeights
    from ._sym import symCheck

    from ._util import setClassArgs, _setDefaultFits, _getFitInds

    # from pemtk.util import hvPlotters  # Imports full module, include setup routines.

    def __init__(self, matE = None, data = None, ADM = None, backend = 'afblmXprod', **kwargs):

        self.__notebook__ = isnotebook()
        self.setClassArgs(locals())  # Set attribs from passed args

        super().__init__(**kwargs)  # Run rest of __init__ from base class - MAY BE PROBLEMATIC WITH ARB **kwargs here, should set separate dict for this?

        # Set dict to hold selection options. May want to propagate back to base class.
        # For using as kwargs set these by function? Might be a bit of a headache, is there a better way?
        # May also want to split by dataType?
        # NOTE: see self.setSubset() for usage, but should propagate elsewhere too.
        self.selOpts = {'matE':{},
                        'slices':{}}

        # 22/08/22 Added fitOpts mainly for backend selection, but should also set other params here (possibly from self.selOpts?)
        self.fitOpts = {'backend':self.backends(backend=backend)}


        # self.subset = {}  # Set dict to hold data subset for fitting
        self.subKey = 'subset'  # Set subKey to use existing .data structure for compatibility with existing functions!
        # self.resultsKey = 'results'  # Set resultsKey to use existing .data structure for compatibility with existing functions!
        self.data[self.subKey] = {}  # Set empty to ensure setMatEFit will always run, even on blank init.

        # Init fit params to None
        # Currently set to push current fit to base structure, and keep all fits in .data[N] dicts.
        self.lmmu = None
        self.basis = None
        self.fitData = None
        self.fitInd = 0

        # Set dict for file handling
        # Currently just a log - appends to list on file IO.
        self.files = {'filesIn':[], 'filesOut':[], 'dataPaths':[], 'filesOutFailed':[]}

        # Quick hack for optional Holoviews import - may not be best way to do this?
        # Issue without binding to self is scope of sub-module/imported methods (e.g. from _analysis.py).
        # Methods within this file ARE within the import's scope however.
        if hvPlotters.hvFlag:
            self.hv = hvPlotters.hv
            self.setPlotDefaults = hvPlotters.setPlotDefaults
            self.setPlotDefaults()
        else:
            self.hv = False
            self.setPlotDefaults = None

        # Same for Seaborn.
        if snsFlag:
            self.sns = sns
        else:
            self.sns = False


    # def setFitSubset(self, thres = None, selDims = None, thresDims = None, sq = True, drop = True):
    #     """
    #     Threshold and subselect on matrix elements.
    #
    #     Wrapper for :py:func:`epsproc.matEleSelector`.
    #
    #     """
    #
    #     self.subset = matEleSelector(self.dsMatE, thres = thres, inds = selDims, dims = thresDims, sq = sq, drop = drop)

    def backends(self, backend = None):
        """
        Set backends & select for use.

        Pass backend = 'name of backend' to select.
        Pass None to return dict of available backends.
        If a function is passed it is returned directly.

        """

        # Define allowed backend fns. by cagegory
        backDict = {'af':{'name':'afblm','afblmXprod':afblmXprod, 'keyDim':'t'},
                    'mf':{'name':'mfblm','mfblmXprod':mfblmXprod, 'keyDim':'pol'}}

        # Set selected backend if passed.
        if backend is not None:
            if isinstance(backend, str):
                if self.verbose['sub']:
                    print(f"* Setting fitting backend to {backend}")
                return backDict[backend[0:2]][backend]

            # If func is passed, just return directly.
            if callable(backend):
                return backend

        else:
            return backDict


# ************* Set & select data

    def setData(self, keyExpt = None, keyData = None):
        """
        Data prototype - this will be used to set experimental data to the master structure.
        For now, set data by passing, and computational data can be set for testing, by passing a key.
        This basically assumes that the expt. provides AFBLMs.

        TO CONSIDER:
        - Data format, file IO for HDF5?
        - Routines to read VMI images and process etc, planned for experimental code-base.
        - Further simulation options, e.g. add noise etc., for testing.

        """

        if isinstance(keyData, str):
            self.data[keyExpt] = {'AFBLM':self.data[keyData]['AFBLM']}
        else:
            self.data[keyExpt] = {'AFBLM':keyData}


    def setWeights(self, wConfig = None, keyExpt = None, keyData = None, **kwargs):
        """
        Wrapper for setting weights for/from data. Basically follows self.setData, with some additional options.

        Will set self.data[keyExpt]['weights'] from existing data if keyData is a string, or from keyData as passed otherwise.

        Parameters
        ----------

        wConfig : optional, str, default = None
            Additional handling for weights.
            - 'poission', set Poissionian weights to match data dims using self.setPoissWeights()
            - 'errors', set weights as 1/(self.data[keyExpt]['weights']**2)


        """

        # Basic setup, as per setData()
        if isinstance(keyData, str):
            self.data[keyExpt]['weights'] = self.data[keyData]['weights']
        else:
            self.data[keyExpt]['weights'] = keyData

        # Additional data handling & options...
        if wConfig == 'poission':
            self.data[keyExpt]['weights'] = setPoissWeights(1.0, self.data[keyExpt]['AFBLM'].shape)
        elif wConfig == 'errors':
            self.data[keyExpt]['weights'] = 1/(self.data[keyExpt]['weights']**2)






    def setADMs(self, **kwargs):
        """Thin wrapper for ep.setADMs(), pass args & set returns to self.data['ADM']"""

        ADMX = setADMs(**kwargs)

        # self.ADM = {'ADMX':ADMX, **kwargs}
        self.data['ADM'] = {'ADM':ADMX}  # Set in main data structure


    def setPolGeoms(self, **kwargs):
        """Thin wrapper for ep.setPolGeoms(), pass args & set returns to self.data['pol']"""

        pol = setPolGeoms(**kwargs)

        # self.ADM = {'ADMX':ADMX, **kwargs}
        self.data['pol'] = {'pol':pol.swap_dims({'Euler':'Labels'})}  # Set in main data structure, force to Labels as dim.
        # self.data['pol'] = pol


#     def setSubset(self, thres = None, selDims = None, thresDims = None, sq = True, drop = True):
    def setSubset(self, dataKey, dataType, sliceParams = None, subKey = None, resetSelectors = False, **kwargs):
        """
        Threshold and subselect on matrix elements.

        Wrapper for :py:func:`epsproc.Esubset` and :py:func:`epsproc.matEleSelector`, to handle data array slicing & subselection from params dict.

        Subselected elements are set to self.data[subKey][dataType], where subKey defaults to self.subKey (uses existing .data structure for compatibility with existing functions!)

        To additionally slice data, set dict of parameters sliceParams = {'sliceDim':[start, stop, step]}

        To reset existing parameters, pass resetSelectors = True.

        To do: better slice handling - will likely have issues with different slices for different dataTypes in current form.

        """
        if subKey is None:
            subKey = self.subKey
        else:
            self.subKey = subKey

        # Arg checks
        if (dataType not in self.selOpts.keys()) or resetSelectors:
            self.selOpts[dataType] = {} # Set empty if missing or to reset.

        for key, value in kwargs.items():
            # self.selOpts['matEleSelector'][key] = value
            self.selOpts[dataType][key] = value  # Set any values if passed


        # Init data['subset']
        if subKey not in self.data.keys():
            self.data[subKey] = {}

        self.data[subKey][dataType] = self.data[dataKey][dataType]


        # Slice on multiple dims if specified, wrap Esubset for this
        # Note stacked dims may cause issues... see http://xarray.pydata.org/en/stable/indexing.html#multi-level-indexing
        #   e.g. data.data['orb5']['matE'].unstack('LM').sel(**{'l':slice(1,3)}) is OK, but data.data['orb5']['matE'].sel(**{'l':slice(1,3)}) is not!
        #        data.data['orb5']['matE'].sel(LM={'l':slice(1,3), 'm':0}) is OK
        # Similar results for loc:
        #        data.data['orb5']['matE'].loc[{'l':slice(1,3), 'm':0}]   is OK, keeps LM
        #        data.data['orb5']['matE'].loc[{'l':slice(1,3)}]   drops l coord
        # if sliceParams is None:
        #     sliceParams = self.d
        #
        if sliceParams is not None:
            for key, value in sliceParams.items():
                self.selOpts['slices'][key] = value
                self.data[subKey][dataType] = self.Esubset(key = subKey, dataType = dataType, Erange = value, Etype = key)

#         self.subset = matEleSelector(self.data[key]['matE'], thres = thres, inds = selDims, dims = thresDims, sq = sq, drop = drop)

        self.data[subKey][dataType] = matEleSelector(self.data[subKey][dataType], **self.selOpts[dataType])


        if self.verbose:
            print(f"Subselected from dataset '{dataKey}', dataType '{dataType}': {self.data[subKey][dataType].size} from {self.data[dataKey][dataType].size} points ({self.data[subKey][dataType].size/self.data[dataKey][dataType].size * 100 :.2f}%)")



# *********** Methods to setup fit
    def setMatEFit(self, matE = None, paramsCons = 'auto', refPhase = 0, colDim = 'it', verbose = 1):
        """
        Convert an input Xarray into (mag,phase) array of matrix elements for fitting routine.

        Parameters
        ----------

        matE : Xarray
            Input set of matrix elements, used to set allowed (l,m,mu) and input parameters.

        paramsCons : dict, optional, default = 'auto'
            Input dictionary of constraints (expressions) to be set for the parameters.
            See https://lmfit.github.io/lmfit-py/constraints.html
            If 'auto', parameters will be set via self.symCheck()

        refPhase : int or string, default = 0
            Set reference phase by integer index or name (string).
            If set to None (or other types) no reference phase will be set.

        colDims : dict, default = 'it'
            Quick hack to allow for restacking via ep.multiDimXrToPD, this will set to cols = 'it', then restack to 1D dataframe.
            This should always work for setting matE > fit parameters, but can be overridden if required.

            This is convienient for converting to Pandas > lmfit inputs, but should redo directly from Xarray for more robust treatment.
            For ePS matrix elements the default should always work, although will drop degenerate cases (it>1). but shouldn't matter here.
            TODO:
               - make this better, support for multiple selectors.
               - For eps case, matE.pd may already be set?


        Returns
        -------
        params : lmfit parameters object
            Set of fitting parameters.

        lmmu : dict
            List of states and mappings from states to fitting parameters (names & indexes).


        29/06/21: Adapted to use 'it' on restack, then set to single-column with dummy dim. No selection methods, use self.setSubset() for this.

        """
        # Quick and ugly wrap args for class - should tidy up here!
        if matE is None:
            # matE = self.subset
            matE = self.data[self.subKey]['matE']

        params = Parameters()
        lmmu = {}  # Dict to hold param names & reindexing for abs/phase to complex mappings.

        # Try passing dict, https://lmfit.github.io/lmfit-py/parameters.html#lmfit.parameter.Parameters.add_many
        # ACTUALLY, not quite right - need tuples per parameter, not dict
        # params.add_many(paramDict)

        # Try list of tuples: having issues with type?
        # params.add_many(testTlist)

        # NOPE - doesn't like tuples as name
        # params.add_many(((0, 1, 0, 0), -2.3170597),
        #                  ((0, 3, 0, 0), 1.1057219))

        # NOPE - invalid name
        # params.add_many(('(0, 1, 0, 0)', -2.3170597),
        #                  ('(0, 3, 0, 0)', 1.1057219))

        # NOPE - complex not supported (presumably), gives type error: TypeError: '>' not supported between instances of 'complex' and 'float'
        # params.add_many(('one', -2.317060+1.358736j),
        #                  ('two', 1.1057219+1.358736j))


        # OK
        # params.add_many(('one', -2.3170597),
        #                  ('two', 1.1057219))


        # Try as list of params with sensible naming... OK
    #     testMatE = data.data['subset']['matE'] #.drop('Sym')  # Subselect before conversion?
        # testMatE = data.data['subset']['matE'].sel({'Eke':1.1},drop=False) #.drop('Sym')  # Setting single Eke here is annoying, as it messes up PD conversion - no way to keep singleton dim???

        # # Using PD conversion routine works, although may have issues with singleton dims again - should set suitable dummy dim here?
        # # pdTest, _ = ep.multiDimXrToPD(testMatE, colDims='Eke', dropna=True)
        # pdTest, _ = multiDimXrToPD(matE, colDims=colDim, dropna=True, squeeze = False)
        # # pdTest, _ = ep.multiDimXrToPD(testMatE, colDims='Sym', dropna=True, squeeze = False)
        #
        # pdTest = pd.DataFrame(pdTest.stack(colDim))  # Stack to 1D format and force to DF
        pdTest = self.pdConvSetFit(matE = matE, colDim = colDim)  # Functional version.

        # Select column from pd dataset - NOW ASSUMED ABOVE, and use .flatten() below to force to 1 column/dim.
    #     col = 1.1
        # pdSub = pdTest[Eke].dropna()
        # pdSub = pdTest.dropna()
        # dataNP = pdSub.to_numpy()  # Drop Eke with col label - probably not best way to do it here.
        # lmmuList = pdSub.index  #.to_numpy  # This returns list of tuples

        dataNP = pdTest.to_numpy().flatten()
        lmmuList = pdTest.index  # This returns list of tuples

        lmmu['Index'] = lmmuList
        nMatE = lmmu['Index'].size

        if verbose:
            print(f"Set {nMatE} complex matrix elements to {2*nMatE} fitting params, see self.params for details.")

        # Convert list of labels to short str format
        # Parameter names [a-z_][a-z0-9_]*, so replace - sign with
    #     lmmuListStr = ['_'.join(str(ele).replace('-','n') for ele in sub) for sub in lmmuList]
        lmmuListStr = lmmuListStrReformat(lmmuList)
        lmmu['Str'] = lmmuListStr

        # Generate (l,m) labels
        lmList = lmmu['Index'].droplevel(['Cont','Targ','Total','mu','it'])  # NOTE this assumes dims, should change to slice?
        strLabels = [','.join(str(ele) for ele in sub) for sub in lmList]

        # Generate map {full names : lm lables}
        lmmu['lmMap'] = dict(zip(lmmu['Str'], strLabels))

    #     return lmmu

        # Set mapping of labels (n, item, abs(item), phase(item)) to use for param remapping
        lmmu['map'] = [(item, [n, n+nMatE], 'm_' + item, 'p_' + item) for n, item in enumerate(lmmu['Str'])]

    #     # Set params with numbered elements
    #     [params.add(f"mag_{n}", np.abs(item), min = 1e-4, max = 5.0) for n, item in enumerate(dataNP)]   # Set min == theshold? Not sure if this will cause issues later.
    #     [params.add(f"phase_{n}", np.angle(item), min=-np.pi, max=np.pi) for n, item in enumerate(dataNP)]

    #     # Set ref phase
    #     params['phase_0'].vary = False

        # Set params with named elements
    #     [params.add(f"m_{lmmuListStr[n]}", np.abs(item), min = 1e-4, max = 5.0) for n, item in enumerate(dataNP)]   # Set min == theshold? Not sure if this will cause issues later.
    #     [params.add(f"p_{lmmuListStr[n]}", np.angle(item), min=-np.pi, max=np.pi) for n, item in enumerate(dataNP)]

        # Set parameters with mapping
        [params.add(lmmu['map'][n][2], np.abs(item), min = 1e-4, max = 5.0) for n, item in enumerate(dataNP)]   # Set min == theshold? Not sure if this will cause issues later.
        [params.add(lmmu['map'][n][3], np.angle(item), min=-np.pi, max=np.pi) for n, item in enumerate(dataNP)]

        # Set any passed constraints
        if paramsCons is not None:
            if paramsCons == 'auto':
                paramsCons, _ = self.symCheck(pdTest = pdTest, colDim = colDim)

                if verbose:
                    print("Auto-setting parameters.")


            for k,v in paramsCons.items():
                if k in params.keys():
                    params[k].set(expr = v)
                else:
                    print(f'*** Warning, parameter {k} not present, skipping constraint {k} = {v}')

        # Set ref phase
    #     if refPhase is not None:
        if isinstance(refPhase,int):
            params[lmmu['map'][0][3]].vary = False
        elif isinstance(refPhase,str):
            params[refPhase].vary = False


    #     return params, lmmuList
        # return params, lmmu
        self.params = params
        self.lmmu = lmmu

        # Set also to data dict for easy reference later (handy for quick pickle save/load)
        self.data[self.subKey]['params'] = self.params
        self.data[self.subKey]['lmmu'] = self.lmmu

        if self.verbose:
            if self.__notebook__:
                display(self.params)
            else:
                print(self.params)



    def reconParams(self, params = None, lmmuList = None):
        """
        Convert parameters object > Xarray for tensor multiplication.

        VERY UGLY!  Should be a neater way to do this with existing Xarray and just replace/link values (i.e. pointers)?

        ... but it does work.

        """

        if params is None:
            params = self.params
        if lmmuList is None:
            lmmuList = self.lmmu

        # Return params to complex vals...
        paramList = []
        v = params.valuesdict()

        # For numbered list
    #     for n in range(0,lmmuList.shape[0]):
    #         paramList.append(v[f'mag_{n}']*np.exp(v[f'phase_{n}']*1j))

        # For named list
    #     for n in range(0,lmmuList.shape[0]):
        # List with reformatting
    #     paramList = [v[f'm_{item}']*np.exp(v[f'p_{item}']*1j) for item in lmmuListStrReformat(lmmuList)]
        # Dict mapping form
        paramList = [v[item[2]]*np.exp(v[item[3]]*1j) for item in lmmuList['map']]


        # xr.DataArray(paramList) # OK
        # xr.DataArray(paramList, coords = lmmuList)  # Not OK

    #     pdRecon = pd.DataFrame(paramList, lmmuList)  # OK

        pdRecon = pd.DataFrame(paramList, lmmuList['Index'])  # Dict form

        xrRecon = xr.DataArray(pdRecon)

        return xrRecon.unstack().stack({'LM':['l','m'], 'Sym':['Cont','Targ','Total']}).squeeze(drop=True)


    def randomizeParams(self):
        """Set random values for self.params."""

        for item in self.params:
        #     item.value = np.random.random()
            if self.params[item].vary:
                self.params[item].value = np.random.random()


# ************ MAIN FITTING ROUTINES

    # def afblmMatEfit(self, matE = None, lmmuList = None, data = None, basis = None, ADM = None, pol = None, selDims = {}, thres = None, thresDims = 'Eke', lmModelFlag = False, XSflag = True, **kwargs):
    def afblmMatEfit(self, matE = None, data = None, lmmuList = None, basis = None, ADM = None, pol = None, resetBasis = False,
                        selDims = {}, thres = None, thresDims = 'Eke', lmModelFlag = False, XSflag = True,
                        weights = None, backend = None, debug = False, **kwargs):
        """
        Wrap :py:func:`epsproc.geomFunc.afblmXprod` for use with lmfit fitting routines.

        Parameters
        ----------

        matE : Xarray or lmfit Parameters object
            Matrix elements to use in calculation.
            For Parameters object, also require lmmuList to convert to Xarray.
            If not passed, use self.data[self.subKey]['matE'].

        data : Xarray, optional, default = None
            Data for fitting.
            If set, return residual.
            If not set, return model result.

        lmmuList : list, optional, default = None
            Mapping for paramters.
            Uses self.lmmu if not passed.

        basis : dict, optional, default = None
            Pre-computed basis set to use for calculations.
            If not set try to use self.basis, or passed set of ADMs.
            NOTE: currently defaults to self.basis if it exists, pass resetBasis=True to force overwrite.

        ADM : Xarray
            Set of ADMs (alignment parameters). Not required if basis is set.

        pol : Xarray
            NOTE: currently NOT used for :py:func:`epsproc.geomFunc.afblmXprod`
            Set of polarization geometries (Euler angles). Not required if basis is set.
            (If not set, defaults to ep.setPolGeoms())

        resetBasis : bool, optional, default=False
            Force self.basis overwrite with updated values.
            NOT YET IMPLEMENTED

        weights : int, Xarray or np.array, optional, default = None
            Weights to use for residual calculation.
            - If set, return np.sqrt(weights) * residual. (Must match size of data along key dimension(s).)
            - If None, try to use use self.data[self.subKey]['weights'].
              If that is not found, or is None, an unweighted residual will be returned.

            For bootstrap sampling, setting Poissonian weights can be used, see https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Poisson_bootstrap
            Use self.setWeights() for this, e.g. weights = rng.poisson(weights, data.t.size)
            To use uncertainties from the data, set weights = 1/(sigma^2)


        backend : function, optional, default = None
            UPDATED 21/08/22 - now default = None, uses self.fitOpts['backend']
            Set at class init, see also self.backends().

            Testing 12/08/22
            Supports backend = afblmXprod or mfblmXprod, and test OK with latter.
            NOTE - when passing fn externally, it may need to be defined in base namespace. (Default case uses locally defined function.)
            E.g.
                data.afblmMatEfit(backend = ep.mfblmXprod) should be OK following `import epsproc as ep`
                data.afblmMatEfit(backend = mfblmXprod) will fail



        NOTE:

        - some assumptions here, will probably need to run once to setup (with ADMs), then fit using basis returned.
        - Currently fitting abs matrix elements and renorm Betas. This sort-of works, but gives big errors on |matE|. Should add options to renorm matE for this case, and options for known B00 values.

        TODO:

        - Consolidate weights to main data structure.
            04/05/22: added to basis return as basis['weights'], may want to pipe back to self.data[self.subKey]['weights'], or just set elsewhere?
        - More sophisticated bootstrapping methods, maybe with https://github.com/smartass101/xr-random and https://arch.readthedocs.io/en/latest/index.html

        21/08/22: now with improved backend handling, working for AF and MF case.
        12/08/22: testing for MF fitting. Initial tests for case where BASIS PASSED ONLY, otherwise still runs AF calc.
        02/05/22: added weights options and updated docs.

        """

        if backend is None:
            backend = self.fitOpts['backend']

        if debug:
            print(f"Running fits with {backend.__name__}")

        # Quick and ugly wrap args for class - should tidy up here!
        if matE is None:
            # matE = self.subset
            matE = self.data[self.subKey]['matE']

        if lmmuList is None:
            lmmuList = self.lmmu

        # if data is None:
        #     # data = self.data
        #     data = self.fitData

        # if resetBasis:
        #     basis =

        if basis is None:
            # if hasattr(self,'basis') and (not resetBasis) and (self.basis is not None):
            if (not resetBasis) and (self.basis is not None):
                basis = self.basis

            # TODO: set other optional params here
            else:
                if (ADM is None) and ('ADM' in self.data[self.subKey].keys()):
                    # ADM = self.ADM['subset']
                    ADM = self.data[self.subKey]['ADM']

                if (pol is None) and ('pol' in self.data[self.subKey].keys()):
                    pol = self.data[self.subKey]['pol']

        if isinstance(matE, Parameters):
            # Set matE from passed params object
            matE = self.reconParams(matE, lmmuList)

        # Set to use XS in fit (default), or renorm betas only (B00=1)
        # Note this is only used for no basis case, then will be set in basis dictionary.
        if XSflag:
            BLMRenorm = 0
        else:
            BLMRenorm = 1

        # Run AFBLM calculation; set basis if not already set.
        if basis is None:
            BetaNormX, basis = backend(matE, AKQS=ADM,   # FOR AF ONLY
                                           RX=pol,  # FOR MF ONLY - RX removed in ePSproc v1.3.0 for AF - not required/valid for AF calcs.
                                           thres = thres, selDims = selDims, thresDims=thresDims,
                                           basisReturn = 'ProductBasis', BLMRenorm = BLMRenorm, **kwargs)

    #         return BetaNormX, basis
            # if resetBasis:
            #     self.basis = basis

        else:
            # Pass **basis here to allow for passing generically through fitting routine and easy/flexible unpacking into afblmXprod()
            BetaNormX = backend(matE, **basis,
                                               thres = thres, selDims = selDims, thresDims=thresDims, basisReturn = 'BLM', **kwargs)

    #         return BetaNormX

        # Setup weights if required - v1 with flexibility
        # if weights is not None:
        #     # Poissonian weights
        #     if isinstance(weights, int) or isinstance(weights, float):
        #         # rng = np.random.default_rng()
        #         # weights = rng.poisson(weights, BetaNormX.shape)
        #         weights = setWeights(weights, BetaNormX.shape)
        #
        #     # Use existing settings if True
        #     elif weights == True:
        #         weights = self.data[self.subKey]['weights']
        #
        # # Add to basis for return if required.
        # basis['weights'] = weights

        # Weights v2 - just set as per pol, ADM etc. from class.
        # If None will be skipped later in any case.
        if weights is None:
            if 'weights' in self.data[self.subKey].keys():
                weights = self.data[self.subKey]['weights']

        if data is not None:
            if weights is None:
                return (np.abs(BetaNormX - data)).values.flatten()   # Need to return NP array of residual for lmfit minimize fitting routine

            # Poissonian weights
            # elif isinstance(weights, int) or isinstance(weights, float):
            #     rng = np.random.default_rng()
            #     weights = rng.poisson(weights, BetaNormX.shape)
                #
                # return (np.sqrt(weights) * np.abs(BetaNormX - data)).values.flatten()

            # Weights as passed
            else:
                return (np.sqrt(weights) * np.abs(BetaNormX - data)).values.flatten()

        else:
            if lmModelFlag:
                return BetaNormX.values.flatten()   # Need to return NP array of model output for lmfit Model routines

            else:
                return BetaNormX, basis


    # Wrap fit routine
    def fit(self, fcn_args = None, fcn_kws = None, fitInd = None, keepSubset = False, **kwargs):
        """

        Wrapper to run lmfit.Minimizer, for details see https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.Minimizer

        Uses preset self.params for parameters, and self.data[self.subKey] for data.

        Default case runs a Levenberg-Marquardt minimization (method='leastsq'), using :py:func:`scipy.optimize.least_squares`, see the
        `Scipy docs for more options <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html>`_, using
        the AF fitting model :py:func:`epsproc.geomFunc.afblmXprod` calculation routine. For MF fitting backend set fcn_kws['backend'] = ep.geomFunc.mfblmXprod


        Parameters
        -----------
        fcn_args : tuple, optional, default = None
            Positional arguments to pass to the fitting function.
            If None, will be set as (self.data[self.subKey]['AFBLM'], self.lmmu, self.basis)

        fcn_kws : dict, optional, default = {}
            Keyword arguments to pass to the fitting function.
            For MF fitting backend set fcn_kws['backend'] = ep.geomFunc.mfblmXprod


        fitInd : int, optional, default = None
            If None, will use self.fitInd
            For parallel usage, supply explicit fitInd instead of using class var to ensure unique key per fit.

        keepSubset : bool, optional, default = False
            If True, keep a copy of self.data[self.subKey] in self.data[fitInd][self.subKey]

        **kwargs
            Passed to the fitting functions, for options see:

            - For lmfit options and defaults see https://lmfit.github.io/lmfit-py/fitting.html
            - For scipy (lmfit backend) see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html


        18/08/22: debugged for MF fitting case, now can pass MF backend via fcn_kws['backend'] = ep.geomFunc.mfblmXprod

        02/05/22: added **kwags for backends.

        07/09/21: updating for parallel use.
                    Note that main outputs (self.reults etc.) are now dropped. May want to set to last result?

        """
        if fitInd is None:
            fitInd = self.fitInd
            self.fitInd += 1

        if fcn_args is None:
            fcn_args = (self.data[self.subKey]['AFBLM'], self.lmmu, self.basis)

        # 21/08/22 - now set defaults in self.fitOpts
        #            TODO: should handle updates here too, but only used for backend functionality at the moment.
        if fcn_kws is None:
            fcn_kws = self.fitOpts


        # Setup fit
    #     self.minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data['subset']['AFBLM']))  # , fcn_args=(lmmuList, BetaNormX, basis))  # In this case fn. args present as self.arg, but leave data as passed arg.
    #     minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data['subset']['AFBLM'].sel(Eke=1.1), self.lmmu, self.basis))  # Above not working with or without 'self', try explicit args instead... THIS IS WORKING, almost, but not quite passing correct things...

        # minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data['subset']['AFBLM'].sel(Eke=1.1), self.lmmu, self.basis))  # Now working OK, just need to sort data pass/setting.
        minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=fcn_args, fcn_kws=fcn_kws, **kwargs)  # Now working OK, just need to sort data pass/setting.

        # Setup with function version to check arg passing OK - NOTE ORDERING is currently different!
    #     minner = Minimizer(afblmMatEfit, self.params, fcn_args=(self.lmmu, self.data['subset']['AFBLM'], self.basis))

        # Run fit
        # self.result = minner.minimize()
        result = minner.minimize()

        # Check final result
        # BetaNormX, _ = self.afblmMatEfit(matE = result.params, **fcn_kws)
        BetaNormX, _ = self.afblmMatEfit(result.params, None, *fcn_args[1:], **fcn_kws)  # DON'T pass data for BLM return.
                                                                                         # TODO: may want to use all kwargs here for clarity.
        # self.betaFit = BetaNormX
        # self.residual = self.afblmMatEfit(matE = self.result.params, data = self.data[self.subKey]['AFBLM'])
        # residual = self.afblmMatEfit(matE = result.params, data = self.data[self.subKey]['AFBLM'], **fcn_kws)
        residual = self.afblmMatEfit(result.params, *fcn_args, **fcn_kws)

        #************ Push results to main data structure
        # May want to keep multiple sets here?
        # if not (self.resultsKey in self.data.keys()):
        #     self.data[self.resultsKey] = {}
        #     self.fitInd = 0

        # Version with simple numerically-indexed results.
        self.data[fitInd] = {}
        self.data[fitInd]['AFBLM'] = BetaNormX.copy()
        # self.data[fitInd]['residual'] = self.residual.copy()
        # self.data[fitInd]['results'] = copy.deepcopy(self.result)  # Full object copy here.
        self.data[fitInd]['residual'] = residual.copy()
        self.data[fitInd]['results'] = result  #.copy()  # Full object copy here.

        # Keep subset data? Useful in some cases
        if keepSubset:
            self.data[fitInd][self.subKey] = self.data[self.subKey].copy()

        # Add some metadata
        timeString = dt.now().strftime('%Y-%m-%d_%H-%M-%S')

        try:
            self.data[fitInd]['AFBLM'].attrs['jobLabel'] = f"Fit #{fitInd}, ({self.data[fitInd]['AFBLM']['t'].size} t, {self.data[fitInd]['AFBLM']['Labels'].size} pol) points, $\chi^2$={self.data[fitInd]['results'].chisqr}\n {timeString}"
        except:
            self.data[fitInd]['AFBLM'].attrs['jobLabel'] = f"Fit #{fitInd}, ({self.data[fitInd]['AFBLM']['Labels'].size} pol) points, $\chi^2$={self.data[fitInd]['results'].chisqr}\n {timeString}"

        # self.fitInd += 1
