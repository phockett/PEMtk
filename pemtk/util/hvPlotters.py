"""
Holoviews util routines for PEMtk

14/09/21    Quick implementation for fit analysis routines.
            May want to consolidate with epsproc.plot.hvPlotters (basically same idea, but different details)

"""


#*** Plot setup
# Pulled code from tmoDataBase.py
# May already be implemented in some form in ePSproc.

# HV imports
try:
    import holoviews as hv
    from holoviews import opts
    hv.extension('bokeh', 'matplotlib')
    hvFlag = True

except ImportError:
    print("*** Holoviews not found: interactive plots not available (hv backend).")
    hvFlag = False


# Set some default plot options
def setPlotDefaults(fSize = [800,400], imgSize = 500):
    """Basic plot defaults"""

    if hvFlag:
        opts.defaults(opts.Scatter(width=fSize[0], height=fSize[1], tools=['hover'], show_grid=True),
                      opts.Curve(width=fSize[0], height=fSize[1], tools=['hover'], show_grid=True),
                      opts.Image(width=imgSize, frame_width=imgSize, aspect='square', tools=['hover'], colorbar=True),   # Force square format for images (suitable for VMI)
                      opts.HeatMap(width=imgSize, frame_width=imgSize, aspect='square', tools=['hover'], colorbar=True),
                      opts.HexTiles(width=fSize[0], height=fSize[1], tools=['hover'], colorbar=True))
