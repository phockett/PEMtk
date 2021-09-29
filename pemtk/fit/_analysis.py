# PEMtk fitting analysis routines
#
# 13/09/21  v1  Basic codes from dev notebook.
#
# Initial dev code: see https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_analysis_demo_150621-tidy.html
#

import xarray as xr
import numpy as np

from epsproc import matEleSelector, multiDimXrToPD

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

    AFpdLong = AFpd.reset_index().melt(id_vars=['Fit','l','m'])  # This works for pushing to full long-format with col by t


    #*** Set to self
    # loc = locals()
    # self.fitMetrics = {i: loc[i] for i in ('dfLong', 'dfRef', 'AFxr', 'AFpd', 'AFxrRS', 'AFpdLong')}  # White list
    # self.fitMetrics = {i: loc[i] if i not in ('AFstack') for i in loc}  # Black list
    # self.fitMetrics = {i: loc[i] for i in loc}  # All
    # self.fitMetrics = locals()  # All - may need .copy()?
    self.data['fits'] = locals()  # All - may need .copy()?
    del self.data['fits']['self']  # Remove 'self'


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

    key : str, optional, default = 'fit'
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

    # Set plot data
    pData = self.data[key][dataDict][dataType]

    if thres is not None:
        self.thresFits(thres = thres, dataType = dataType, key = key, dataDict = dataDict)

    if ('mask' in self.data[key].keys()) and mask:
        if len(pData) == len(self.data[key]['mask'][dataType]):
            pData = pData[self.data[key]['mask'][dataType]]

            print(f"Mask selected {self.data[key]['mask'][dataType].sum()} results (from {self.data[key]['mask'][dataType].count()}).")

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
