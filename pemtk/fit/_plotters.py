# Wrappers for some ePSproc plotting routines


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
