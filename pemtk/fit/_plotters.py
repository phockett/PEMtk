# Wrappers for some ePSproc plotting routines


def BLMfitPlot(self, keys = ['subset'], dataType='AFBLM', Etype='t', thres=1e-2, col=None, **kwargs):
    """Wrap BLMplot for data + simulation results with default params."""

    # Add current fit results
    if self.fitInd-1 >= 0:
        keys.append(self.fitInd-1)

    # Check keys OK
    keys = self._keysCheck(keys)

    # super().BLMplot(keys=keys, dataType=dataType, thres = thres, Etype=Etype)  # In this case renamed function to avoid messing up other routines using BLMplot
    self.BLMplot(keys=keys, dataType=dataType, thres = thres, Etype=Etype, **kwargs)
