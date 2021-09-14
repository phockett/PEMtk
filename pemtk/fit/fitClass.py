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
from epsproc.geomFunc import afblmXprod

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

class pemtkFit(dataClass):
    """
    Class prototype for pemtkFit class. Dev version builds on dataClass, and adds some basic subselection & fitting capabilities.
    """

    from ._conv import pdConv, pdConvSetFit
    from ._util import setClassArgs
    from ._parallel import multiFit
    from ._plotters import BLMfitPlot, lmPlotFit
    from ._sym import symCheck

    def __init__(self, matE = None, data = None, ADM = None, **kwargs):

        self.__notebook__ = isnotebook()
        self.setClassArgs(locals())  # Set attribs from passed args

        super().__init__(**kwargs)  # Run rest of __init__ from base class - MAY BE PROBLEMATIC WITH ARB **kwargs here, should set separate dict for this?

        # Set dict to hold selection options. May want to propagate back to base class.
        # For using as kwargs set these by function? Might be a bit of a headache, is there a better way?
        # May also want to split by dataType?
        self.selOpts = {'matE':{},
                        'slices':{}}

        # self.subset = {}  # Set dict to hold data subset for fitting
        self.subKey = 'subset'  # Set subKey to use existing .data structure for compatibility with existing functions!
        # self.resultsKey = 'results'  # Set resultsKey to use existing .data structure for compatibility with existing functions!

        # Init fit params to None
        # Currently set to push current fit to base structure, and keep all fits in .data[N] dicts.
        self.lmmu = None
        self.basis = None
        self.fitData = None
        self.fitInd = 0


    # def setFitSubset(self, thres = None, selDims = None, thresDims = None, sq = True, drop = True):
    #     """
    #     Threshold and subselect on matrix elements.
    #
    #     Wrapper for :py:func:`epsproc.matEleSelector`.
    #
    #     """
    #
    #     self.subset = matEleSelector(self.dsMatE, thres = thres, inds = selDims, dims = thresDims, sq = sq, drop = drop)

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
            if paramsCons is 'auto':
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
                        selDims = {}, thres = None, thresDims = 'Eke', lmModelFlag = False, XSflag = True, **kwargs):
        """
        Wrap :py:func:`epsproc.geomFunc.afblmXprod` for use with lmfit fitting routines.

        Parameters
        ----------

        matE : Xarray or lmfit Parameters object
            Matrix elements to use in calculation.
            For Parameters object, also require lmmuList to convert to Xarray.
            If not passed, use self.subset.


        ADM : Xarray.
            Set of ADMs (alignment parameters). Not required if basis is set.

        pol : Xarray
            Set of polarization geometries (Euler angles). Not required if basis is set.
            (If not set, defaults to ep.setPolGeoms())


        basis : optional, default = None
            Use to pass pre-computed basis set.
            NOTE: currently defaults to self.basis if it exists, pass resetBasis=True to force overwrite.

        NOTE:

        - some assumptions here, will probably need to run once to setup (with ADMs), then fit using basis returned.
        - Currently fitting abs matrix elements and renorm Betas. This sort-of works, but gives big errors on |matE|. Should add options to renorm matE for this case, and options for known B00 values.

        """

        # Quick and ugly wrap args for class - should tidy up here!
        if matE is None:
            # matE = self.subset
            matE = self.data[self.subKey]['matE']

        if lmmuList is None:
            lmmuList = self.lmmu

        # if data is None:
        #     # data = self.data
        #     data = self.fitData

        if basis is None:
            # if hasattr(self,'basis') and (not resetBasis) and (self.basis is not None):
            if (not resetBasis) and (self.basis is not None):
                basis = self.basis

            # TODO: set other optional params here
            else:
                if ADM is None:
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
            BetaNormX, basis = afblmXprod(matE, AKQS=ADM,   # RX=pol,  # RX removed in ePSproc v1.3.0 - not required/valid for AF calcs.
                                           thres = thres, selDims = selDims, thresDims=thresDims,
                                           basisReturn = 'ProductBasis', BLMRenorm = BLMRenorm, **kwargs)

    #         return BetaNormX, basis

        else:
            # Pass **basis here to allow for passing generically through fitting routine and easy/flexible unpacking into afblmXprod()
            BetaNormX = afblmXprod(matE, **basis,
                                               thres = thres, selDims = selDims, thresDims=thresDims, basisReturn = 'BLM', **kwargs)

    #         return BetaNormX

        if data is not None:
            return (np.abs(BetaNormX - data)).values.flatten()   # Need to return NP array of residual for lmfit minimize fitting routine

        else:
            if lmModelFlag:
                return BetaNormX.values.flatten()   # Need to return NP array of model output for lmfit Model routines

            else:
                return BetaNormX, basis


    # Wrap fit routine
    def fit(self, fitInd = None):
        """
        Wrapper to run lmfit.Minimizer.

        Uses preset self.* parameters.

        For parallel usage, supply explicit fitInd instead of using class var.

        07/09/21: updating for parallel use.
                    Note that main outputs (self.reults etc.) are now dropped. May want to set to last result?

        """
        if fitInd is None:
            fitInd = self.fitInd
            self.fitInd += 1


        # Setup fit
    #     self.minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data['subset']['AFBLM']))  # , fcn_args=(lmmuList, BetaNormX, basis))  # In this case fn. args present as self.arg, but leave data as passed arg.
    #     minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data['subset']['AFBLM'].sel(Eke=1.1), self.lmmu, self.basis))  # Above not working with or without 'self', try explicit args instead... THIS IS WORKING, almost, but not quite passing correct things...

        # minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data['subset']['AFBLM'].sel(Eke=1.1), self.lmmu, self.basis))  # Now working OK, just need to sort data pass/setting.
        minner = Minimizer(self.afblmMatEfit, self.params, fcn_args=(self.data[self.subKey]['AFBLM'], self.lmmu, self.basis))  # Now working OK, just need to sort data pass/setting.

        # Setup with function version to check arg passing OK - NOTE ORDERING is currently different!
    #     minner = Minimizer(afblmMatEfit, self.params, fcn_args=(self.lmmu, self.data['subset']['AFBLM'], self.basis))

        # Run fit
        # self.result = minner.minimize()
        result = minner.minimize()

        # Check final result
        BetaNormX, _ = self.afblmMatEfit(matE = result.params)
        # self.betaFit = BetaNormX
        # self.residual = self.afblmMatEfit(matE = self.result.params, data = self.data[self.subKey]['AFBLM'])
        residual = self.afblmMatEfit(matE = result.params, data = self.data[self.subKey]['AFBLM'])

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


        # Add some metadata
        timeString = dt.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.data[fitInd]['AFBLM'].attrs['jobLabel'] = f"Fit #{fitInd}, ({self.data[fitInd]['AFBLM']['t'].size} t, {self.data[fitInd]['AFBLM']['Labels'].size} pol) points, $\chi^2$={self.data[fitInd]['results'].chisqr}\n {timeString}"

        # self.fitInd += 1
