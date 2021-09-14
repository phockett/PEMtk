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

def analyseFits(self):
    """
    Collate fit data from multiple runs.

    See https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_analysis_demo_150621-tidy.html for dev code.

    """

    #*** Reformat data from class.

    # Convert fit results to Pandas
    dfLong, dfRef = self.pdConv()

    # Stack AFBLM fit outputs to Xarrays
    # Don't recall the best way to do this - here pass to list then stack directly

    # AFstack = xr.DataArray()
    AFstack = []

    for n in range(0, self.fitInd):
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
    self.fitMetrics = locals()  # All - may need .copy()?


def fitHist(self, bins = None, dataType = 'redchi', pRange = None, backend = 'hv'):
    """
    Basic histogram plot of batch fit results.

    TODO:

    - see TMO-DEV (https://github.com/phockett/tmo-dev) for some generalised plotting methods.
    - Implement, but better, with decorators for data checking & selection.
    - Import chain: currently using util.hvPlotters.py > self.hv for backend, but could also use a decorator here?
    - Data subselection by threshold or range. (Again see TMO-DEV routines for ideas.)

    - Holoviews stuff
       - Fix data subset to plotter, otherwise get full dataset to tooltip.
       - Hist bar options to fix.

    """

    # Set plot data
    pData = self.fitMetrics['dfLong'][dataType]

    # Pandas/Matplotlib histogram
    if backend in ['pd', 'mpl']:
        pData.hist(bins=bins)

    # Holoviews scatter + histograms
    elif backend is 'hv':
        hvObj = self.hv.Scatter(pData.reset_index(), kdims=dataType).hist(dimension=[dataType,'Fit'])

        # Code from showPlot()
        if self.__notebook__:     # and (not returnImg):
            display(hvObj)  # If notebook, use display to push plot.
