# PEMtk fitting analysis routines
#
# 26/11/21  v2  Most functionality now implemented.
#               Very messy, and not well-integrated.
#               TODO:   better integration, esp. options etc. to classes or class dicts.
#
#
# 13/09/21  v1  Basic codes from dev notebook.
#
# Initial dev code: see https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_analysis_demo_150621-tidy.html
#
# TODO: add metadata to PD DataFrames, currently some routines assumwe name and/or try/fail to allow for different types.
#       Need to differentiate between wide/long parameter sets, and fit function outputs.
#       UPDATE: now added basic dTypes in df.attrs['dType'], but may not propagate.
#       ALSO: not yet very general in many cases, due to current structures.

import pprint
import string

import pandas as pd
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt   # For alternative sns routines


from epsproc import matEleSelector, multiDimXrToPD
from epsproc.util.misc import subselectDims

from ._util import phaseCorrection as phaseCorrFunc
from ._util import addColLevel, renameParams

# #*** Plot setup - NOW MOVED TO .util.hvPlotters.py
# NOTE POSSIBLE SCOPE ISSUES HERE - currently resolved by importing to self.hv, but may not be best practice.

# # Pulled code from tmoDataBase.py
# # May already be implemented in some form in ePSproc.
#
# # HV imports
# try:
#     import holoviews as hv
#     from holoviews import opts
#     hv.extension('bokeh', 'matplotlib')
#
# except ImportError:
#     print("*** Holoviews not found: interactive plots not available (hv backend).")
#
#
# # Set some default plot options
# def setPlotDefaults(fSize = [800,400], imgSize = 500):
#     """Basic plot defaults"""
#     opts.defaults(opts.Scatter(width=fSize[0], height=fSize[1], tools=['hover'], show_grid=True),
#                   opts.Curve(width=fSize[0], height=fSize[1], tools=['hover'], show_grid=True),
#                   opts.Image(width=imgSize, frame_width=imgSize, aspect='square', tools=['hover'], colorbar=True),   # Force square format for images (suitable for VMI)
#                   opts.HeatMap(width=imgSize, frame_width=imgSize, aspect='square', tools=['hover'], colorbar=True),
#                   opts.HexTiles(width=fSize[0], height=fSize[1], tools=['hover'], colorbar=True))


# def scopeTest(self):
#     print("Checking isnotebook()")
#     print(isnotebook())


def analyseFits(self, dataRange = None, batches = None):
    """
    Collate fit data from multiple runs.

    Data from self.
    For individual

    See https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_analysis_demo_150621-tidy.html for dev code.

    """

    dataRange = self._setDefaultFits(dataRange)

    #*** Reformat data from class.

    # Convert fit results to Pandas
    dfLong, dfRef = self.pdConv(dataRange = dataRange, batches = batches)

    # Set per-fit metrics
    dfPF = dfLong.droplevel(['Type','pn']).reset_index().drop_duplicates('Fit').drop(['Param','vary','expr','value','stderr'], axis=1).set_index('Fit')
    dfPF.attrs['dType'] = 'Per Fit metrics'

    # For per-fit analysis, set also a "wide" format table
    # Potential issues in this form with updating, so just compose later as necessary for plotting.
    # dfWide = dfLong.reset_index().pivot_table(columns = 'Param', values = ['value'], index=['Fit','Type'],aggfunc=np.sum)


    # Stack AFBLM fit outputs to Xarrays
    # Don't recall the best way to do this - here pass to list then stack directly

    # AFstack = xr.DataArray()
    AFstack = []

    for n in range(dataRange[0], dataRange[1]):
        try:
            AFstack.append(self.data[n]['AFBLM'].expand_dims({'Fit':[n]})) # NOTE [n] here, otherwise gives length not coord
                                                                       # http://xarray.pydata.org/en/stable/generated/xarray.DataArray.expand_dims.html
        except KeyError as e:
            if self.verbose['sub']:
                print(f"*** Missing fit key {n}")


    AFxr = xr.concat(AFstack,'Fit')
    AFxr.attrs['jobLabel'] = "BLM results (all fits)"  # TODO: put some useful info here, date, time, # fits etc.

    AFpd, AFxrRS = multiDimXrToPD(AFxr.squeeze().pipe(np.abs), colDims=['t'],
                                 thres = 1e-4, fillna=True)
    AFpd.attrs['dType'] = 'AF results Wide'
    AFpdLong = AFpd.reset_index().melt(id_vars=['Fit','l','m'])  # This works for pushing to full long-format with col by t
    AFpdLong.attrs['dType'] = 'AF results Long'


    #*** Set to self
    # loc = locals()
    # self.fitMetrics = {i: loc[i] for i in ('dfLong', 'dfRef', 'AFxr', 'AFpd', 'AFxrRS', 'AFpdLong')}  # White list
    # self.fitMetrics = {i: loc[i] if i not in ('AFstack') for i in loc}  # Black list
    # self.fitMetrics = {i: loc[i] for i in loc}  # All
    # self.fitMetrics = locals()  # All - may need .copy()?
    self.data['fits'] = locals()  # All - may need .copy()?
    del self.data['fits']['self']  # Remove 'self'

    self._setWide()   # Set wide format params with current data, will reset later.

    self.fitsReport()  # Generate some default report data

    # Dict for plot objects - TODO, currently set in plotting fncs. Should centralise & set options.
    # self.data[]


def _setData(self, key, dataDict, dataType = None, thres = None, mask = True):
    # Set plot data - full dict
    pData = self.data[key][dataDict]

    # Subset
    if dataType is not None:
        pData = pData[dataType]

    # Threshold
    if thres is not None:
        self.thresFits(thres = thres, dataType = dataType, key = key, dataDict = dataDict)
        # pData = pData[self.data[key]['mask'][dataType]]

    if ('mask' in self.data[key].keys()) and mask:
        try:
            if len(pData) == len(self.data[key]['mask'][dataType]):
                pData = pData[self.data[key]['mask'][dataType]]

                print(f"Mask selected {self.data[key]['mask'][dataType].sum()} results (from {self.data[key]['mask'][dataType].count()}).")
        except KeyError:
            if self.verbose['sub']:
                print(f"Mask not set for dataType = {dataType}.")

    return pData


def _setWide(self, indexDims = ['Fit','Type','chisqrGroup','redchiGroup'], valueDims = ['value'],
                key = 'fits', dataDict = 'dfLong', dataWide = 'dfWide',
                dataIn = None, returnFlag = False):
    """
    Set "wide" format data table from current dfLong array, with optional index dims.

    Defaults to self.data[key][dataDict], or pass dataIn to use this instead and return wide-form data directly.

    NOTE: current form checks indexDims for valid subselection. May want to add smarter dim matching for groups?

    """
    if self.verbose['sub']:
        print(f"Setting wide-form data self.[{key}][{dataWide}] from self.[{key}][{dataDict}] (as pivot table).")

    if dataIn is not None:
        dataDict = dataIn

    elif self.data[key][dataDict].attrs['dType'].endswith('Long'):
        dataDict = self.data[key][dataDict]


    # if self.data[key][dataDict].attrs['dType'].endswith('Long'):
    #
    #     iDims = subselectDims(self.data[key][dataDict], refDims = indexDims)  # Check dims
    #     vDims = subselectDims(self.data[key][dataDict], refDims = valueDims)  # Check dims
    #
    #     self.data[key][dataWide] = self.data[key][dataDict].reset_index().pivot_table(columns = 'Param', values = vDims, index=iDims,aggfunc=np.sum)
    #     self.data[key][dataWide].attrs['dType'] = 'Params Wide'

    else:
        print(f"Please set or pass long-form DataFrame, self.data[{key}][{dataDict}] has type {self.data[key][dataDict].attrs['dType']}.")
        return 0


    # Ugly but working
    # Note set as str if single-membered list - otherwise always get multindex on output dataframe
    iDims = subselectDims(dataDict, refDims = indexDims)  # Check dims
    # iDims = [item for item in iDims if item in indexDim]  # Force ordering to match input.
    iDims=iDims[0] if len(iDims)==1 else iDims

    vDims = subselectDims(dataDict, refDims = valueDims)  # Check dims
    vDims=vDims[0] if len(vDims)==1 else vDims

    if self.verbose['sub']:
        # print(f"Setting wide-form data self.[{key}][{dataWide}] from self.[{key}][{dataDict}] (as pivot table).")
        print(f"Index(es) = {iDims}, cols = {vDims}")

    dfWide = dataDict.reset_index().pivot_table(columns = 'Param', values = vDims, index=iDims, aggfunc=np.sum)
    dfWide.attrs['dType'] = 'Params Wide'

    if returnFlag:
        return dfWide
    else:
        self.data[key][dataWide] = dfWide


def fitsReport(self, key = 'fits', dataDict = 'dfPF', thres = None, mask = True):
    """
    Generate fit report/metrics, defaults to self.data['fits']['dfPF'].
    Results are printed if self.verbose, and also set to self.fitsSummary.

    Parameters
    ----------
    key : str, optional, default = 'fits'
        Key into main self.data dictionary.

    dataDict : str, optional, default = 'dfPF'
        Dataset to use, from self.data[key].
        Default case is per-fit metrics.

    thres : float, optional, default = None
        Set threshold for subselection, for range [0, thres].
        This is passed to self.thresFits() and sets self.data[key]['mask'].

    mask : bool, optional, default = True
        Use self.data[key]['mask'] to subselect data if set.

    """

    # Set data subset, functional form - may already have better function elsewhere...?
    pData = self._setData(key, dataDict,  thres = thres, mask = mask)

    # Get metrics
    fitReport = {'Fits': pData.index.unique(level='Fit').shape[0],  # More general          # pData.shape[0],  # Basic
             'Success': pData['success'].sum(),   # Only valid for dataType = dfPF at the moment, should generalise.
             'Minima': {'chisqr':pData['chisqr'].min(), 'redchi':pData['redchi'].min()}
            }

    if self.verbose['main']:
        pprint.pprint(fitReport, indent = 2)

    self.fitsSummary = fitReport


# Test final param determination/report by group
def paramsReport(self, key = 'fits', dataDict = 'dfWide', inds = {}, aggList = ['min', 'mean', 'median', 'max', 'std', 'var']):
    """
    Generate parameter report/metrics, defaults to self.data['fits']['dfWide'].
    Results are printed if self.verbose, and also set to self.paramsSummary.

    Parameters
    ----------
    key : str, optional, default = 'fits'
        Key into main self.data dictionary.

    dataDict : str, optional, default = 'dfWide'
        Dataset to use, from self.data[key].
        Default case is per-fit metrics.

    inds : dict, optional, default = {}
        Set of indexs to subselect from, as dictionary items.
        E.g. xs = {'redchiGroup':'C'} will select group 'C'.

    aggList : list, optional, default = ['min', 'mean', 'median', 'max', 'std', 'var']
        List of aggregator functions to use.
        These are passed to Pandas.agg(), a list of common functions can be found at https://pandas.pydata.org/docs/user_guide/basics.html#descriptive-statistics

    TODO: consolidate indexing methods with other functions & extend to thesholds and cross-ref (column) values.

    """


    # subset = self.data[key][dataDict]

    # Basic subselection (indexes) - needs dim checking
    # for k,v in inds.items():
    #     subset = subset.xs(v, level = k)

    subset = self._setData(key, dataDict)
    subset = self._subsetFromXS(selectors = inds, data = subset)

    # Stats per dataType with groups
    self.paramsSummary = {}
    self.paramsSummary['data'] = subset
    self.paramsSummary['count'] = subset.count()

    # Default decribe()
    self.paramsSummary['desc'] = subset.groupby('Type').describe().T

    # Custom agg
    self.paramsSummary['agg'] = subset.groupby('Type').agg(aggList).T
    self.paramsSummary['agg'].index.rename(['Param','Agg'], inplace=True)

    if self.verbose['main']:
        print("Set parameter stats to self.paramsSummary.")
        if self.__notebook__:
#             display(self.paramsReport['desc'])
            display(self.paramsSummary['agg'])




def paramsCompare(self, params = None, ref = None, phaseCorr = True,
                    phaseCorrParams = {}):
    """
    Compare extracted parameter set with reference data.

    NOTE: currently assumes self.paramsSummary and self.params for aggregate fit results & reference data.

    Parameters
    ----------
    params : pd.DataFrame, optional, default = None
        Fit parameters to tabulate.
        Will use self.paramsSummary in default case (and run self.paramsReport() if missing).

    ref : pd.DataFrame, optional, default = None
        Reference parameter set to compare with.
        Will use self.data['fits']['dfRef'] in default case (and attempt to set this if missing).

    phaseCorr : bool, optional, default = True
        Run phase correction routine for reference parameters.

    phaseCorrParams : dict, optional, default = {}
        Pass dictionary to additionally set parameters for phaseCorrection() method.
        Default cases runs with {'dataDict':'dfRef', 'dataOut':'dfRefPC', 'useRef':False}, these parameters will update the defaults.
        Note - these params are only used if phaseCorr = True.


    TODO:

    - Better dim handling.
    - Generalize to compare any pair of parameter sets. (Just need to loop over param sets and check attrs labels.)

    """

    if params is None:

        if not hasattr(self, 'paramsSummary'):
            print(f"Missing self.paramsSummary, running self.paramsReport() to generate with defaults.")
            self.paramsReport()

        params = self.paramsSummary.copy()


    if ref is None:

#         if not hasattr(self, 'params'):
#             self.setMatEFit(paramsCons = {})


#         dfRef = self.params.copy()

        if ('dfRef' in self.data['fits'].keys()) and (self.data['fits']['dfRef'] is not None):
            dfRef = self.data['fits']['dfRef']  # use this if already set?
        else:
            dfRef = self.pdConvRef()  # Set ref from self.params, will also be set if missing.
            self.data['fits']['dfRef'] = dfRef

    else:
        dfRef = ref

    if phaseCorr:
        # Set defaults & update with any passed args
        phaseCorrArgs = {'dataDict':'dfRef', 'dataOut':'dfRefPC', 'useRef':False}
        phaseCorrArgs.update(phaseCorrParams)

        # Run phaseCorrection()
        self.phaseCorrection(**phaseCorrArgs)
        dfRef = self.data['fits']['dfRefPC']

    # For wide-form data as returned by phaseCorrection()
    if dfRef.attrs['dType'].endswith('Wide'):
        dfRefRestack = dfRef.T  #.reset_index()  #.pivot(columns = 'Type')  #.reset_index()
        dfRefRestack.columns = dfRefRestack.columns.reorder_levels(['Type','Fit'])

    # For long-form data as set in original params table.
    else:
        # UGLY reshape - should be in existing routines already somewhere...? NOTE ISSUE WITH PIVOT > MULTIINDEX COLS not required here.
        # dfRefRestack = dfRef.reset_index().set_index(['Param']).drop(['Fit','pn','stderr', 'vary', 'expr'], axis=1).pivot(columns = 'Type')
        dfRefRestack = dfRef.reset_index().set_index(['Param'])[['Type','value']].pivot(columns = 'Type')  # Cleaner to just keep specified cols
        dfRefRestack = dfRefRestack.droplevel(None, axis=1)  # Drop extra level, defaults to None
        addColLevel(dfRefRestack)
        dfRefRestack.columns = dfRefRestack.columns.reorder_levels(['Type','Dataset'])
        # addColLevel(dfRefRestack, newCol = 'input',  names = ['Dataset','Type','Agg'])

    dataFit = params['agg'].unstack().copy()  #.xs('mean', level='Agg')  #.unstack()  # Without xs to keep other Agg values.
    # addColLevel(dataFit, newCol = 'Fit')

    # **** Sort & append.
    # TODO: tidy this up, lots of .T since this was all developed quasi-independently with slightly different data formats.
    dataMerge = dataFit.merge(dfRefRestack, on = 'Param').T

    # print(dataMerge)
    # Add level for values vs. %
    # See also ._util.addColLevel()
    newLevel = 'dType'
    baseName = 'num'  # TODO: set data type here based on inputs? Real/imag/abs/float?
    pcName = '%'
    dataMerge[newLevel] = baseName
    dataMerge.set_index(newLevel, append=True, inplace=True)
    dataMerge = dataMerge.T
    # print(dataMerge)

    # # Set differences - no trivial way to do this?
    # for item in dataMerge.columns.get_level_values('Type').unique():
    #     try:
    #         dataMerge[item,'diff'] = dataMerge[item]['mean'] - dataMerge[item]['ref']   #.diff(-1)  #['m']
    #     except:
    #         pass   # Crude error-handling for missing dims

    # Update with additional values
    for item in dataMerge.columns.get_level_values('Type').unique():
        try:
            dataMerge[item,'diff',baseName] = dataMerge[item]['mean'] - dataMerge[item]['ref']   #.diff(-1)  #['m']
            dataMerge[item,'diff/std', pcName] = np.abs(dataMerge[item,'diff']/dataMerge[item]['std'] *100)
            # dataMerge[item,'diff/std %']['Cat'] = ['%']

            for subitem in ['std','diff']:
                dataMerge[item, f'{subitem}', pcName] = np.abs(dataMerge[item][subitem]/dataMerge[item]['mean'] * 100)   #.diff(-1)  #['m']

        except:
            pass   # Crude error-handling for missing dims

    # dfT = dataMerge.T.sort_index(level = 'Type')  #.reindex(['mean','ref','diff','std'], level = 'Agg')  # NOTE reindex requires all levels here.
    # dfT


    # Final sorting and reformat
    # print(dataMerge)
    dfOut = dataMerge.T.sort_index(level = 'Type').reindex(['mean','ref','diff','std','diff/std'], level = 'Agg')
    # dfOut = dataMerge.T.sort_index(level = 'Type').reindex(['mean','ref','diff','std'], level = 'Agg')
    # print(dfOut)
    dfOut.index.set_names(['Type','Source',newLevel], inplace = True)

    # Set outputs
    self.paramsSummaryComp = dfOut

    if self.verbose['main']:
        print("Set parameter comparison to self.paramsSummaryComp.")
        if self.__notebook__:
#             display(self.paramsReport['desc'])
            display(self.paramsSummaryComp)


#************* Classifications & transformations

def classifyFits(self, key = 'fits', dataDict = 'dfPF', dataType = 'redchi', group = None, bins = None, labels = None,
                    plotHist = True, propagate = True):
    """
    Classify fit result sets (DataFrame) based on chisqr or redchi values.

    Parameters
    ----------

    bins : int or list, default = None
        Bins setting for classifier
        - Set as None for default case, will bin by (min - min*0.05, min*5, 10)
        - Set as int to define a specific number of (equally spaced) bins, for (min - min*0.05, min*5, numbins)
        - Set as a list [start,stop] or [start,stop,bins] for specific range.
        - Set as list (>3 elements) to define specific bin intervals.

    dataType : str, default = 'redchi'
        DataType to classify.
        (Currently only supports a single dataType.)

    key : str, optional, default = 'fits'
        Key into main self.data dictionary.

    dataDict : str, optional, default = 'dfPF'
        Dataset to use, from self.data[key].
        Default case is per-fit metrics.

    group : str, optional, default = None
        Name for classifier column in dataframe.
        If None, defaults to {dataType}Group

    labels : list, optional, default = None
        Specify names for group memebers.
        Defaults to alphabetic labels.

    plotHist : bool, optional, default = True
        Plot histogram of results & show tabular summary (from self.data[key][group])

    propagate : bool, optional, default = True
        propagate classifications to other data types?

    """

    # Set data subset, functional form - may already have better function elsewhere...?
    pData = self._setData(key, dataDict, dataType = dataType,  thres = None, mask = False)

    min = self.fitsSummary['Minima'][dataType]

    # Set binning - case for None or number of bins only
    defaultBins = 10
    if (bins is None) or isinstance(bins,int):
        bins = np.linspace(min - min * 0.05, min * 5, (bins if isinstance(bins,int) else defaultBins))
    elif len(bins) == 2:
        bins = np.linspace(bins[0], bins[1], defaultBins)
    elif len(bins) == 3:
        bins = np.linspace(bins[0], bins[1], bins[2])

    if labels is None:
        labels = list(string.ascii_uppercase[0:len(bins)-1])

    if group is None:
        group = dataType + 'Group'

    self.data[key][dataDict][group] = pd.cut(pData, bins = bins, labels = labels)
    # self.data[key][group] = pd.DataFrame(np.array([bins[0:-1], bins[1:]]).T, index = labels, columns=['Min','Max'])  # Set array of bins
    self.data[key][group] = self.data[key][dataDict].groupby(group).describe()  # Set descriptive array + bins
    self.data[key][group]['Min'] = bins[0:-1]
    self.data[key][group]['Max'] = bins[1:]

    self.data[key][group].attrs['dType'] = 'Group metrics'
    self.data[key][group].attrs['bins'] = bins
    self.data[key][group].attrs['group'] = group

    # Check result - note this is ordered by first appearance in the dataFrame, unless sorted first. Also NaNs are dropped.
    # data.data['fits']['dfPF']['pGroups'].sort_values().dropna().hist(bins=10)
    if plotHist:
        # self.fitHist(key = key, dataDict = dataDict, dataType = group, backend = 'pd')  # This code isn't quite general enough as yet!
                                                                                          # HV part fails with categorical data at the moment.
        self.data[key][dataDict][group].sort_values().dropna().hist()

        if self.__notebook__:
            display(self.data[key][group])

    # TODO: tidy this up! All special cases at the moment.
    #       Esp. dfLong, which currently reruns pd.cut without error checks. Ugh. Horrible.
    #       BETTER APPROACH: just set self.data[key][group] as index or mask into other dataframes? (Or similar - don't need to set in all existing DFs anyway.)
    if propagate:
        for table in self.data[key].keys():
            if table not in [dataDict, group]:
                try:
                    if isinstance(self.data[key][table], pd.DataFrame):
                        if self.data[key][table].attrs['dType'].endswith('Long'):
                            # Specific update for wide params table.
                            if self.data[key][table].attrs['dType'] == 'Params Long':
                                self.data[key][table][group] = pd.cut(self.data[key][table][dataType], bins = bins, labels = labels)  # Ugh, just do it again for Params Long, otherwise multiindex gets dropped.
                                self._setWide(key = key, dataDict = table, indexDims = ['Fit','Type',group])

                            else:
                                # Try propagating... OK, although odd results for dfLong & dfWide (due to multiindex?)
                                attrs = self.data[key][table].attrs
                                # self.data[key][table] = self.data[key][table].merge(self.data[key][dataDict][group].reset_index(), on='Fit', how='left')
                                self.data[key][table] = self.data[key][table].merge(self.data[key][dataDict][group], on='Fit', how='left')

                                # Keep attrs - merge seems to remove these sometimes otherwise
                                self.data[key][table].attrs = attrs

                            if self.verbose['sub']:
                                print(f"Set {group} for data frame {table}.")

                # Generic exception handling
                except Exception as e:
                    if self.verbose['sub']:
                        print(f"Couldn't set {group} for data frame {table}. Error {type(e)}: {e.args}.")


def phaseCorrection(self, key = 'fits', dataDict = 'dfLong', dataOut = 'dfWide', dataRef = 'dfRef', useRef = True, returnFlag = False, **kwargs):
    """
    Wrapper for ._util.phaseCorrection() (functional form).

    Parameters
    ----------

    key : str, default = 'fits'
        Data key for analysis dataframe.

    dataDict : str, default = 'dfLong'
        Data dict for analysis dataframe.
        Note default case uses `self.data[key][dataDict]`

    dataOut : str, default = 'dfWide'
        Output dict key for phase-corrected dataframe.

    dataRef : str, default = 'deRef'
        Reference dict key for phase.
        Default case uses `self.data[key][dataRef]`
        Note this is ONLY USED IF useRef = True is set.

    useRef : bool, default = True
        Use reference phase from self.data[key][dataRef]?
        Otherwise ref phase will be set to phasesIn.columns[0] (as per :py:func:`pemtk.fit._util.phaseCorrection`)

    returnFlag : bool, default = True
        If True return phase-corrected data.
        If False, set data to self.data[key][dataOut]

    **kwargs
        Passed to :py:func:`pemtk.fit._util.phaseCorrection`


    NOTE: this currently only sets phaseCorrected data in wide-form dataset, self.data[key][dataOut]. May want to push to long-form too? (Otherwise this will be lost by self._setWide().)

    TODO: tidy up options here, a bit knotty.

    """

    # Work with copy and set phase corr data to this...
    dataIn = self.data[key][dataDict].copy()
    # dataIn.index = dataIn.index.set_levels(['m','pc'], level = 'Type')  # Set 'pc' Type
    dataInWide = self._setWide(dataIn = dataIn, returnFlag = True)  # Set to wide form
    refDims = dataInWide.index.names

    if useRef:
        dfRef = self.data[key][dataRef]
    else:
        dfRef = None

    dfOut = phaseCorrFunc(dataInWide, dfRef = dfRef, **kwargs)  # Run phase corr routine

    # Update data with new values as 'pc' Type - CASE FOR LONG DATA
    # dfOut.index = dfOut.index.set_levels(['m','pc'], level = 'Type')  # Set 'pc' Type
    # # dfOut.name = 'value'
    # dfOut = pd.concat([dfOut, dataInWide]).sort_index()  # This sort-of works, get all values but duplicate m
    # dfOut = dfOut[~dfOut.index.duplicated(keep='first')]  # Remove any duplicated indexers

    # For wide data - this will be phase only in current case (as returned from phaseCorrFunc())
    # Add Type & reindex (lost in XC in phaseCorrFunc, but may change in future)
    dfOut['Type'] = 'pc'
    # dfOut = dfOut.reset_index().set_index(['Fit','Type'])  # OK, returns new DF
    dfOut = dfOut.reset_index().set_index(refDims)  # Match to original df
    dfOut = pd.concat([dfOut, dataInWide]).sort_index()  # NOTE - this seems to mess up with Multiindex IF dim ordering is different. UGH. HORRIBLE.
                                            # Update: Dim ordering should now be enforced in self._setWide() for dataInWide.
    dfOut.attrs['dType'] = 'Params Wide'

    if self.__notebook__:
        display(dfOut)

    if returnFlag:
        return dfOut
    else:
        self.data[key][dataOut] = dfOut





# ***************************************************************************************
#
# Plotters for fit data/analysis.
#
#

def fitHist(self, bins = 'auto', dataType = 'redchi', key = 'fits', dataDict = 'dfPF',
            thres = None, mask = True, binRange = None, backend = 'hv', plotDict = 'plots'):
    """
    Basic histogram plot of batch fit results.

    Parameters
    ----------
    bins : str, int or list, default = 'auto'
        Bins setting for histogram, essentially as per Numpy routine https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html
        - Set as string for various auto options.
        - Set as int to define a specific number of (equally spaced) bins.
        - Set as list to define specific bin intervals.
        NOTE: some combinations currently not working with 'hv' backend.

    dataType : str, default = 'redchi'
        DataType to histogram.
        (Currently only supports a single dataType.)

    key : str, optional, default = 'fits'
        Key into main self.data dictionary.

    dataDict : str, optional, default = 'dfPF'
        Dataset to use, from self.data[key].
        Default case is per-fit metrics.

    thres : float, optional, default = None
        Set threshold for plotting, for range [0, thres].
        This is passed to self.thresFits() and sets self.data[key]['mask'].
        For more control use binRange setting.

    mask : bool, optional, default = True
        Use self.data[key]['mask'] to subselect data if set.

    binRange : list, optional, default = None
        Specify range for binning.
        Note this is only used by HV plotter, and will override `bins` settings for auto types.
        Specify `bins = int` and `binRange = [start, stop]` for full control.

    backend : str, optional, default = 'hv'
        Specify backend:
        - 'hv' for Holoviews
        - 'pd' or 'mpl' for Pandas.hist()

    plotDict : str, optional, default = 'plots'
        For hv case, return plot object & data to self.data[plotDict] as ['fitHistPlot'] and ['fitHistData']

    Notes:

    - Data to plot is specified by self.data[key][dataDict][dataType].
    - Threshold value sets mask, this will overwrite existing selection mask.
    - If self.data[key]['mask'] exists, this will be used if mask = True.

    TODO:

    - see TMO-DEV (https://github.com/phockett/tmo-dev) for some generalised plotting methods.
    - Implement, but better, with decorators for data checking & selection.
    - Import chain: currently using util.hvPlotters.py > self.hv for backend, but could also use a decorator here?
    - Data subselection by threshold or range. (Again see TMO-DEV routines for ideas.)
    - Fix binning issues with certain cases.

    - Holoviews stuff
       - Fix data subset to plotter, otherwise get full dataset to tooltip.
       - Hist bar options to fix. UPDATE: now set to bins='auto' as default, which works well.
       - See hv.help(histogram) or http://holoviews.org/user_guide/Transforming_Elements.html for more.

    """

    # # Set plot data
    # pData = self.data[key][dataDict][dataType]
    #
    # if thres is not None:
    #     self.thresFits(thres = thres, dataType = dataType, key = key, dataDict = dataDict)
    #
    # if ('mask' in self.data[key].keys()) and mask:
    #     if len(pData) == len(self.data[key]['mask'][dataType]):
    #         pData = pData[self.data[key]['mask'][dataType]]
    #
    #         print(f"Mask selected {self.data[key]['mask'][dataType].sum()} results (from {self.data[key]['mask'][dataType].count()}).")

    # Functional form - may already have better function elsewhere...?
    pData = self._setData(key, dataDict, dataType = dataType,  thres = thres, mask = mask)

    # Clean up data to per-fit properties (unless these are required)
    # TODO: check dims instead of assuming here!
    if (dataType not in ['Type','pn']) and (dataDict is not 'dfPF'):
        pData = pData.drop_duplicates().droplevel(['Type','pn'])

    # Pandas/Matplotlib histogram
    if backend in ['pd', 'mpl']:
        pData.hist(bins=bins)

    # Holoviews scatter + histograms
    elif backend is 'hv':
        # Set bins, treat int bins as separate parameter num_bins for HV.
        # For HV: http://holoviews.org/reference_manual/holoviews.operation.html#holoviews.operation.histogram
        binsHV = (bins if not isinstance(bins,int) else None)
        num_bins = (bins if isinstance(bins,int) else None)

        # bin_range = (bins if (isinstance(bins,list) & (len(bins)==2)) else None) # Use independent pRange parameter for this
        if binRange is not None:
            binsHV = None  # bin_range only applies if bins = None.
        # else:

        # Create plot object.
        # print(bins, binsHV, num_bins, binRange)
        hvObj = self.hv.Scatter(pData.reset_index(), kdims=dataType).hist(dimension=[dataType,'Fit'], bins = binsHV, num_bins = num_bins, bin_range = binRange)

        # Set output
        if not plotDict in self.data.keys():
            self.data[plotDict] = {}

        self.data[plotDict]['fitHistData'] = pData
        self.data[plotDict]['fitHistPlot'] = hvObj

        # Code from showPlot()
        if self.__notebook__:     # and (not returnImg):
            display(hvObj)  # If notebook, use display to push plot.

def _mergePFLong(self, pData, key, dataPF, hue, hRound):
    """
    Merge per-fit data to existing dataset & force to long form.

    Used for Seaborn plotting routines with hue mapping.
    """

    hDims = subselectDims(pData, refDims = hue)
    if hDims:
        # Annoying special case.
        if hue != 'Fit':
            pDataLong = pData.reset_index().melt(id_vars=['Fit',hue])
        else:
            pDataLong = pData.reset_index().melt(id_vars=['Fit'])

    # If hDims is missing, assume it's in daPF - reformat & merge.
    else:
        hData = self.data[key][dataPF][hue]  # TODO: add dim checks here.

        pDataLong = pData.reset_index().melt(id_vars=['Fit'])
        pDataLong = pDataLong.merge(hData, on = ['Fit'], how='left')

        if hRound is not None:
            # Round/rebin hue data to specified dp.
            pDataLong[hue] = pDataLong[hue].apply(np.round, decimals = hRound)

    return pDataLong


def corrPlot(self, key = 'fits', dataDict = 'dfWide', hue = 'redchiGroup', hRound = None,
            dataType = None, level = None, sel = None, selLevel = 'redchiGroup',
            dataPF = 'dfPF', plotDict = 'plots', remap = None,
            backend = 'sns', pairgrid = False, **kwargs):
    """
    Similar to paramPlot(), but set for correlation matrix plotter.

    This requires wide-form parameters data self.data['fits']['dfWide'].
    Two levels of selection are currently supported (index data only, NOT columns)

    **kwargs are passed to Seaborn's pairplot routine, https://seaborn.pydata.org/generated/seaborn.pairplot.html

    TODO: if numerical data columns are added for hue mapping they may result in additional plots too.

    TODO: add HV gridmatrix + linked brushing: http://holoviews.org/user_guide/Linked_Brushing.html

    TODO: FIX HORRIBLE SELECTION ROUTINES.

    """

    # Set plot data from dict.
    # dataPlot = data.data[key][dataDict].xs(dataType,level='Type')
    pData = self._setData(key, dataDict)  #, dataType = dataType)  #,  thres = thres, mask = mask)

    # Basic single-level selection
    # TODO: update selectors as looped dict specs?
    # SEE _util._subsetFromXS
    if dataType is not None:
        pData = pData.xs(dataType, level=level)

    # Further subselection if specified
    if sel is not None:
        pData = pData.xs(sel, level = selLevel)  # , drop_level=False)

    # Handle hue
    # pDataLong = self._mergePFLong(pData, hue)  # Function currently sets long form, here need to retain wide!
    hDims = subselectDims(pData, refDims = hue)

    if hDims:
        pData = pData.reset_index()
        # Annoying special case.
        if hue != 'Fit':
            pData = pData.drop('Fit', axis=1)  # Use existing column in dataframe

    else:
        pData = pData.reset_index()  #.drop('Fit', axis=1)
        hData = self.data[key][dataPF][hue]
        pData = pData.merge(hData, on = ['Fit'], how='left').drop('Fit', axis=1)  # Add per-fit coloumn.

    if hRound is not None:
        # Round/rebin hue data to specified dp.
        pData[hue] = pData[hue].apply(np.round, decimals = hRound)

    # Rename params for plotting?
    if remap is not None:
        pData = renameParams(pData, self.lmmu[remap])  # Remap columns
        # pDataLong.replace({'Param':self.lmmu[remap]}, inplace=True)


    # Create plot
    # if self.sns:
    #     g = self.sns.pairplot(pData, hue = hue, **kwargs)
    #     # sns.pairplot(data.data['fits']['dfWide'], hue='redchiGroup')
    #     # g.set_xticklabels(rotation=-60)
    #
    #     # Code from showPlot()
    #     if self.__notebook__:     # and (not returnImg):
    #         display(g)  # If notebook, use display to push plot.
    #
    # else:
    #     print("Seaborn not loaded, paramPlot() not available.")

    # Create plot
    if self.sns and (backend == 'sns'):
        if (pData[hue].drop_duplicates().count() > 30) or pairgrid:
            # PairGrid - much faster for continuous data cmap
            # From https://stackoverflow.com/questions/60859082/is-there-a-way-to-color-the-points-in-a-python-seaborn-pairplot-with-a-quantitat
            g = self.sns.PairGrid(pData, hue = hue, **kwargs).map(plt.scatter)

        else:
            g = self.sns.pairplot(pData, hue = hue, **kwargs)
            # g.set_xticklabels(rotation=-60)

    # elif self.hv and (backend == 'hv'):
    #  TODO: implement gridmatrix here. See http://holoviews.org/user_guide/Linked_Brushing.html


    else:
        g = ''
        # print("Please set 'sns' (Seaborn or Holoviews backend loaded, paramPlot() not available.")
        print(f"Plotting backed '{backend}' unavailable.")

    # Code from showPlot()
    if self.__notebook__ and g:     # and (not returnImg):
        display(g)  # If notebook, use display to push plot.

    if not plotDict in self.data.keys():
        self.data[plotDict] = {}

    self.data[plotDict]['corrData'] = pData
    self.data[plotDict]['corrPlot'] = g



def paramPlot(self, dataType = 'm', level = 'Type', sel = None, selLevel = 'redchiGroup',
            hue = None, hRound = 7, x='Param', y='value',
            key = 'fits', dataDict = 'dfWide', dataPF = 'dfPF', plotDict = 'plots',
            # thres = None, mask = True,
            hvType = None, remap = None,
            backend = 'sns', returnFlag = False):
    """
    Basic scatter-plot of parameter values by name/type.

    Currently supports Seaborn for back-end only, and requires wide-form dataDict as input.

    TODO:
    - better and more concise dim handling for multi-level selection. Integrate to single dict of selectors? (See tmo-dev?)
    - Box/violin plot options. Also option to drop scatter plot in these cases (now partially implemented for HV only).
    - HV support?
       - Basic support now in place, but cmapping needs some work for non-cat data. SEE NOTES ELSEWHERE!
       - Also breaks for subselection case unless another hue dim is set.
       - 18/11/21: better, but messy, support now in place. Includes Violin & BoxWhisker options.
       - TODO: implement grouping and/or holomap for extra dims.

    Currently: have `selLevel` and `hue`, which must be different in general.
    - `sel` and `selLevel` define subselection by a value in a column, e.g. sel = 'E', selLevel = 'redchiGroup' for values E in column 'selLevel'
    - `hue` specifies hue mapping for Seaborn plot, which must be a column name.
    - If `hue` is not in input data, it will be taken from the per-fit dataframe.

    Ref: Seaborn catplot, https://seaborn.pydata.org/generated/seaborn.catplot.html
    Ref: HV scatter,

    For usage notes see https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy.html

    """

    # ADDTIONAL NOTES from PKG version testing.
    # TODO: fix issues with hue mapping!
    #
    # data.paramPlot(backend='hv')  # OK
    # data.paramPlot(dataType = 'm', hue = 'redchiGroup', backend='hv')  # OK
    #
    #
    # Nope - need grouping or overlay probably, or accidentally messing up DF here?
    # TypeError: '<' not supported between instances of 'float' and 'str'
    # data.paramPlot(dataType = 'm', hue = 'redchi', backend='hv')
    # data.paramPlot(dataType = None, hue = 'redchiGroup', backend='hv')
    # data.paramPlot(dataType = None, hue = 'redchi', backend='hv')

    # Set plot data from dict.
    # dataPlot = data.data[key][dataDict].xs(dataType,level='Type')
    pData = self._setData(key, dataDict)  #, dataType = dataType)  #,  thres = thres, mask = mask)

    # TODO: update selectors as looped dict specs?
    # SEE _util._subsetFromXS
    # Subselect with XS
    if dataType is not None:
        pData = pData.xs(dataType, level=level)

    # Further subselection if specified
    if sel is not None:
        pData = pData.xs(sel, level = selLevel)  # , drop_level=False)

    # Melt to long form for sns.catplot functionality
    # TODO: check & fix any non-consistent dim handling here!
    if hue is None:
        hue = selLevel

    # hDims = subselectDims(pData, refDims = hue)
    # if hDims:
    #     pDataLong = pData.reset_index().melt(id_vars=['Fit',hue])
    #
    # # If hDims is missing, assume it's in daPF - reformat & merge.
    # else:
    #     hData = self.data[key][dataPF][hue]  # TODO: add dim checks here.
    #
    #     pDataLong = pData.reset_index().melt(id_vars=['Fit'])
    #     pDataLong = pDataLong.merge(hData, on = ['Fit'], how='left')
    #
    #     if hRound is not None:
    #         # Round/rebin hue data to specified dp.
    #         pDataLong[hue] = pDataLong[hue].apply(np.round, decimals = hRound)

    pDataLong = self._mergePFLong(pData, key, dataPF, hue, hRound)  # Functionalised version of above.

    # Output current plot data for ref.
    if not plotDict in self.data.keys():
        self.data[plotDict] = {}

    # Rename params for plotting?
    if remap is not None:
        pDataLong.replace({'Param':self.lmmu[remap]}, inplace=True)

    self.data[plotDict]['paramData'] = pDataLong

    # Create plot
    if self.sns and (backend == 'sns'):
        g = self.sns.catplot(x=x, y=y, hue = hue, data = pDataLong)  # pGroups + scatter plot - this shows groupings better
        g.set_xticklabels(rotation=-60)

    elif self.hv and (backend == 'hv'):
        g = self.hv.Scatter(pDataLong, kdims=x, vdims=[y,hue])   # For hv case need to include hue as vdims (or will be dropped). hv to PD class may be cleaner?
        g.opts(color=hue, size=10, jitter=0.4, alpha=0.5, colorbar = True, cmap='coolwarm', legend_position = 'right')  # Some reasonable options, although should pass as **kwargs.
                                                                                                                        # NOTE: for cat data use legend_position = 'right' to force legend out of plot axes.
                                                                                                                        #       For numerical data the cbar appears correctly.
                                                                                                                        # See http://holoviews.org/_modules/holoviews/plotting/bokeh/element.html

        # Add overlays if specified
        # Not stacking order to foreground scatter plot.
        if hvType == 'box':
            ov = self.hv.BoxWhisker(pDataLong, kdims=x, vdims=y)
            g = ov*g
        elif hvType == 'violin':
            ov = self.hv.Violin(pDataLong, kdims=x, vdims=y)
            ov.opts(inner='stick')
            g = ov*g

    else:
        g = ''
        # print("Please set 'sns' (Seaborn or Holoviews backend loaded, paramPlot() not available.")
        print(f"Plotting backed '{backend}' unavailable.")

    # Code from showPlot()
    if self.__notebook__ and g:     # and (not returnImg):
        display(g)  # If notebook, use display to push plot.

    # Output current plot data for ref.
    # self.data['plotData'] = pDataLong
    # self.data['plot'] = g
    # self.data[plotDict]['paramData'] = pData
    self.data[plotDict]['paramPlot'] = g

    if returnFlag:
        return pDataLong

# def corrPlot()
