# Additional plotting methods for pemtk.fitClass
#

import numpy as np
from epsproc import matEleSelector, plotTypeSelector


#******************* Wrappers for some ePSproc plotting routines


def BLMfitPlot(self, keys = None, dataType='AFBLM', Etype='t', thres=1e-2, col=None, **kwargs):
    """
    Wrap BLMplot for data + simulation results with default params.

    TODO:
    - better plotting (HV?).
    - fix legend & colour mapping.
    """

    if keys is None:
        keys = [self.subKey]  # Default to subset data

        # Add current fit results
        if self.fitInd-1 >= 0:
            keys.append(self.fitInd-1)

    if keys is 'all':
        keys = [self.subKey]
        keys.append(list(range(0, self.fitInd)))  # Plot all fits

    # # Check keys OK - not required since BLMplot() will check this
    # keys = self._keysCheck(keys)
    #
    # print(keys)

    # super().BLMplot(keys=keys, dataType=dataType, thres = thres, Etype=Etype)  # In this case renamed function to avoid messing up other routines using BLMplot
    # self.BLMplot(keys=keys, dataType=dataType, thres = thres, Etype=Etype, col=col, **kwargs)  # Pass all keys - single plot, same parameters

    # Quick hack for marker types
    for key in keys:
        # print(key)
        if key is self.subKey:
            self.BLMplot(keys=key, dataType=dataType, thres = thres, Etype=Etype, col=col, marker = 'x', linestyle = 'dashed', **kwargs)
        else:
            self.BLMplot(keys=key, dataType=dataType, thres = thres, Etype=Etype, col=col, **kwargs)


def lmPlotFit(self, keys = None, dataType='AFBLM', Etype='t', thres=1e-2, **kwargs):
        """
        Wrap lmPlot for data + simulation results with default params.

        """

        if keys is None:
            keys = [self.subKey]  # Default to subset data

            # Add current fit results
            if self.fitInd-1 >= 0:
                keys.append(self.fitInd-1)

        self.lmPlot(keys = keys, dataType=dataType, Etype=Etype, thres=thres, **kwargs)


#*********************** General plotters

def BLMsetPlot(self, key = 'fits', dataDict = 'AFxr', agg = True, ref = True,
                overlay = ['l','m'], pType = 'r',
                thres = 1e-3, sel = None,  xDim = None, sq = True, drop = True,
                unstack = True, plotDict = 'plots'):
    """
    Plot sets of BLM results from Xarray datasets with Holoviews.

    For plotting individual datasets with more control, see BLMfitPlot().

    TODO:
    - add Seaborn plotting options.
    - Streamline, should be able to use recursively to stack additional plots...?

    Parameters
    ----------

    agg : bool, default = True
        If True, define reduced data as hv.reduce(['Fit'], np.mean, spreadfn=np.std)
        TODO: more options here.

    ref : bool, default = True
        If True, include original fitted data in plots.
        TODO: more options here.

    """

    # Reduce data using matEleSelector functionality
    dataRed = matEleSelector(self.data[key][dataDict], thres = thres, inds = sel, dims = xDim, sq = sq, drop=drop)

    # May be required depending on dataType and package versions.
    if unstack:
        dataRed = dataRed.unstack()

    # Set plot data type
    daPlot = plotTypeSelector(dataRed, pType = pType, axisUW = xDim)
    hvDS = self.hv.Dataset(daPlot.rename('BLM'))   # TODO: pull name from dataset?

    # SET HV Dataset & plot objects (or use .hvplot?)
    if agg:
        hvDS = hvDS.reduce(['Fit'], np.mean, spreadfn=np.std)

        # Spread over fits + inputs with marker
        hvPlot = hvDS.to(self.hv.Spread, kdims = xDim).overlay(overlay)

    else:
        # Basic line plots
        hvPlot = hvDS.to(self.hv.Curve, kdims = xDim).overlay(overlay)


    # Set ref data - currently use input data only here
    if ref:
        # dataRef = self.data[key][dataDict]
        # dataRef = self.data[self.subKey]['AFBLM']
        dataRef = matEleSelector(self.data[self.subKey]['AFBLM'], thres = thres, inds = sel, dims = xDim, sq = sq, drop=drop)

        if unstack:
            dataRef = dataRef.unstack()

        dataRef = plotTypeSelector(dataRef, pType = pType, axisUW = xDim)
        refDS = self.hv.Dataset(dataRef.rename('BLM-inputs'))   # TODO: pull name from dataset?

        # Add to plot object
        refPlot = refDS.to(self.hv.Curve, kdims = xDim).opts(line_dash='dashed') * refDS.to(self.hv.Scatter, kdims = xDim).opts(marker='cross', size=15)

        # Ugly... but need to ensure matching dims. May be a neater way to force/check this.
        # if agg:
        #     hvPlot = (hvPlot * refPlot).overlay(overlay)
        # else:
        #     hvPlot = hvPlot * refPlot.overlay(overlay)
        hvPlot = hvPlot * refPlot.overlay(overlay)

    # if overlay:
    #     hvPlot = hvPlot.overlay(overlay)  # This fails on Fits if refDS is set, since it's an unmatched dim. But OK if set independently for both datasets, and overlay([]) just passes.

    # Set output
    if not plotDict in self.data.keys():
        self.data[plotDict] = {}

    self.data[plotDict]['BLMsetData'] = hvDS
    self.data[plotDict]['BLMsetRef'] = refDS
    self.data[plotDict]['BLMsetPlot'] = hvPlot

    # Code from showPlot()
    if self.__notebook__:
        display(hvPlot)  # If notebook, use display to push plot.
