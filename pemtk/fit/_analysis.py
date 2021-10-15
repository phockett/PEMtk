# PEMtk fitting analysis routines
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

from epsproc import matEleSelector, multiDimXrToPD
from epsproc.util.misc import subselectDims

from ._util import phaseCorrection as phaseCorrFunc

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


def analyseFits(self, dataRange = None):
    """
    Collate fit data from multiple runs.

    Data from self.
    For individual

    See https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_analysis_demo_150621-tidy.html for dev code.

    """

    dataRange = self._setDefaultFits(dataRange)

    #*** Reformat data from class.

    # Convert fit results to Pandas
    dfLong, dfRef = self.pdConv(dataRange = dataRange)

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
        AFstack.append(self.data[n]['AFBLM'].expand_dims({'Fit':[n]})) # NOTE [n] here, otherwise gives length not coord
                                                                       # http://xarray.pydata.org/en/stable/generated/xarray.DataArray.expand_dims.html

    AFxr = xr.concat(AFstack,'Fit')
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


def _setData(self, key, dataDict, dataType = None, thres = None, mask = True):
    # Set plot data - full dict
    pData = self.data[key][dataDict]

    # Subset
    if dataType is not None:
        pData = pData[dataType]

    # Threshold
    if thres is not None:
        self.thresFits(thres = thres, dataType = dataType, key = key, dataDict = dataDict)

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


def classifyFits(self, key = 'fits', dataDict = 'dfPF', dataType = 'redchi', group = None, bins = None, labels = None,
                    plotHist = True, propagate = True):
    """
    Classify fit result sets (DataFrame) based on chisqr or redchi values.

    Parameters
    ----------

    bins : int or list, default = None
        Bins setting for classifier
        - Set as None for default case, will bin by (min - min*0.05, min*5, 10)
        - Set as int to define a specific number of (equally spaced) bins, for (min - min*0.05, min*10, numbins)
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
                                self.data[key][table][group] = pd.cut(self.data[key][table][group], bins = bins, labels = labels)  # Ugh, just do it again for Params Long, otherwise multiindex gets dropped.
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
    Wrapper for _util.phaseCorrection() (functional form).
    """

    # Work with copy and set phase corr data to this...
    dataIn = self.data[key][dataDict].copy()
    # dataIn.index = dataIn.index.set_levels(['m','pc'], level = 'Type')  # Set 'pc' Type
    dataInWide = self._setWide(dataIn = dataIn, returnFlag = True)  # Set to wide form

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
    # Add Type & reindex
    dfOut['Type'] = 'pc'
    dfOut = dfOut.reset_index().set_index(['Fit','Type'])  # OK, returns new DF
    dfOut = pd.concat([dfOut, dataInWide]).sort_index()  # NOTE - this seems to mess up with Multiindex IF dim ordering is different. UGH. HORRIBLE.
                                            # Update: Dim ordering should now be enforced in self._setWide() for dataInWide.
    # dfOut.attrs['dType'] = 'Params Wide'

    if self.__notebook__:
        display(dfOut)

    if returnFlag:
        return dfOut
    else:
        self.data[key][dataOut] = dfOut


def fitHist(self, bins = 'auto', dataType = 'redchi', key = 'fits', dataDict = 'dfPF',
            thres = None, mask = True, binRange = None, backend = 'hv'):
    """
    Basic histogram plot of batch fit results.

    Parameters
    ----------
    bins : str, int or list, default = 'auto'
        Bins setting for histogram, essentially as per Numpy routine https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html
        - Set as string for various auto options.
        - Set as int to define a specific number of (equally spaced) bins.
        - Set as list to define specific bin intervals.

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

    Notes:

    - Data to plot is specified by self.data[key][dataDict][dataType].
    - Threshold value sets mask, this will overwrite existing selection mask.
    - If self.data[key]['mask'] exists, this will be used if mask = True.

    TODO:

    - see TMO-DEV (https://github.com/phockett/tmo-dev) for some generalised plotting methods.
    - Implement, but better, with decorators for data checking & selection.
    - Import chain: currently using util.hvPlotters.py > self.hv for backend, but could also use a decorator here?
    - Data subselection by threshold or range. (Again see TMO-DEV routines for ideas.)

    - Holoviews stuff
       - Fix data subset to plotter, otherwise get full dataset to tooltip.
       - Hist bar options to fix. UPDATE: now set to bins='auto' as default, which works well.

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
        bins = (bins if not isinstance(bins,int) else None)
        num_bins = (bins if isinstance(bins,int) else None)

        # bin_range = (bins if (isinstance(bins,list) & (len(bins)==2)) else None) # Use independent pRange parameter for this
        if binRange is not None:
            bins = None  # bin_range only applies if bins = None.

        # Create plot object.
        hvObj = self.hv.Scatter(pData.reset_index(), kdims=dataType).hist(dimension=[dataType,'Fit'], bins = bins, num_bins = num_bins, bin_range = binRange)

        # Code from showPlot()
        if self.__notebook__:     # and (not returnImg):
            display(hvObj)  # If notebook, use display to push plot.
