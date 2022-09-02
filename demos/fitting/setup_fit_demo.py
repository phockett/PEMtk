# Setup demo fitting data for PEMtk testing
# Follows approx first half of demo notebook https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html
# Just removed plotting/output functions here for use in further testing.

print('*** Setting up demo fitting workspace and main `data` class object...')
print('For more details see https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html')
print('To use local source code, pass the parent path to this script at run time, e.g. "setup_fit_demo ~/github"')
print('\n\n* Loading packages...')

# A few standard imports...

# +
import sys
import os
from pathlib import Path
# import numpy as np
# import epsproc as ep
# import xarray as xr

from datetime import datetime as dt
timeString = dt.now()
# -

# And local module imports. This should work either for installed versions (e.g. via `pip install`), or for test code via setting the base path below to point at your local copies.

# +
# For module testing, include path to module here, otherwise use global installation
# print(globals())
# NOTE - this currently doesn't pick up preset modPath case in notebook? Not sure why - something incorrect in scoping here.

# With passed arg
args = sys.argv
localFlag = False

if len(args) > 1:
    modPath = Path(args[1])

    # Append to sys path
    sys.path.append((modPath/'ePSproc').as_posix())
    sys.path.append((modPath/'PEMtk').as_posix())

    localFlag = True


try:
    # ePSproc
    import epsproc as ep

    # Import fitting class
    from pemtk.fit.fitClass import pemtkFit

# Default case
# if not 'modPath' in locals():
except ImportError as e:
    if localFlag:
        print(f"\n*** Couldn't import local packages from root {modPath}.")
    else:
        print(f"\n*** Couldn't import packages, are ePSproc and PEMtk installed? For local copies, pass parent path to this script.")


# Set data path
# Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
epDemoDataPath = Path(ep.__path__[0]).parent/'data'



# +
# Set HTML output style for Xarray in notebooks (optional), may also depend on version of Jupyter notebook or lab, or Xr
# See http://xarray.pydata.org/en/stable/generated/xarray.set_options.html
# if isnotebook():
# xr.set_options(display_style = 'html')
# -

# Set some plot options
ep.plot.hvPlotters.setPlotters()

# ### Set & load parameters
#
# There are a few things that need to be configured for a given case...
#
# - Matrix elements (or (l,m,mu) indicies) to use.
# - Alignment distribution (ADMs).
# - Polarization geometry.
#
# In this demo, a real case will first be simulated with computational values, and then used to test the fitting routines.
#
#
# #### Matrix elements
#
# For fit testing, start with computational values for the matrix elements. These will be used to simulate data, and also to provide a list of parameters to fit later.

# +
# Set for ePSproc test data, available from https://github.com/phockett/ePSproc/tree/master/data
# Note this is set here from ep.__path__, but may not be correct in all cases.


# Multiorb data
dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')

print(f'\n* Loading demo matrix element data from {dataPath}...')

data = pemtkFit(fileBase = dataPath, verbose = 0)

# Read data files
data.scanFiles()
data.jobsSummary()

# #### Alignment distribution moments (ADMs)
#
# The class [wraps ep.setADMs()](https://epsproc.readthedocs.io/en/dev/modules/epsproc.sphCalc.html#epsproc.sphCalc.setADMs). This returns an isotropic distribution by default, or values can be set explicitly from a list. Values are set in `self.data['ADM']`.
#
# Note: if this is not set, the default value will be used, which is likely not very useful for the fit!

# Default case
data.setADMs()

# Load time-dependent ADMs for N2 case
# Adapted from ePSproc_AFBLM_testing_010519_300719.m
from scipy.io import loadmat
ADMdataFile = os.path.join(epDemoDataPath, 'alignment', 'N2_ADM_VM_290816.mat')

print(f'\n\n* Loading demo ADM data from {ADMdataFile}...')

ADMs = loadmat(ADMdataFile)

# Set tOffset for calcs, 3.76ps!!!
# This is because this is 2-pulse case, and will set t=0 to 2nd pulse (and matches defn. in N2 experimental paper)
tOffset = -3.76
ADMs['time'] = ADMs['time'] + tOffset

data.setADMs(ADMs = ADMs['ADM'], t=ADMs['time'].squeeze(), KQSLabels = ADMs['ADMlist'], addS = True)
data.data['ADM']['ADM']


# ### Polarisation geometry/ies
#
# This wraps [ep.setPolGeoms](https://epsproc.readthedocs.io/en/dev/modules/epsproc.sphCalc.html#epsproc.sphCalc.setPolGeoms). This defaults to (x,y,z) polarization geometries. Values are set in `self.data['pol']`.
#
# Note: if this is not set, the default value will be used, which is likely not very useful for the fit!
#

data.setPolGeoms()


# +
# # data.setPolGeoms(eulerAngs = [[0,0,0]], labels = ['z'])
# data.setPolGeoms(eulerAngs = [0,0,0], labels = 'z')
# data.data['pol']['pol']  #.swap_dims({'Euler':'Labels'})


# ### Subselect data
#
# Currently handled in the class by setting `self.selOpts`, this allows for simple reuse of settings as required. Subselected data is set to `self.data['subset'][dataType]`, and is the data the fitting routine will use.

# +
# Settings for type subselection are in selOpts[dataType]

print(f'\n* Subselecting data...')

# E.g. Matrix element sub-selection
data.selOpts['matE'] = {'thres': 0.01, 'inds': {'Type':'L', 'Eke':1.1}}
data.setSubset(dataKey = 'orb5', dataType = 'matE')  # Subselect from 'orb5' dataset, matrix elements


# And for the polarisation geometries...
data.selOpts['pol'] = {'inds': {'Labels': 'z'}}
data.setSubset(dataKey = 'pol', dataType = 'pol')

# And for the ADMs...
data.selOpts['ADM'] = {}   #{'thres': 0.01, 'inds': {'Type':'L', 'Eke':1.1}}
data.setSubset(dataKey = 'ADM', dataType = 'ADM', sliceParams = {'t':[4, 5, 4]})


# ## Compute AF-$\beta_{LM}$ and simulate data
#
# With all the components set, some observables can be calculated. For testing, we'll also use this to simulate an experiemental trace...
#
# Here we'll use `self.afblmMatEfit()`, which is also the main fitting routine, and essentially wraps `epsproc.afblmXprod()` to compute AF-$\beta_{LM}$s (for more details, see the [ePSproc method development docs](https://epsproc.readthedocs.io/en/dev/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html)).
#
# If called without reference data, the method returns computed AF-$\beta_{LM}$s based on the input subsets already created, and also a set of (product) basis functions generated - these can be examined to get a feel for the sensitivity of the geometric part of the problem, and will also be used in fitting to limit repetitive computation.

# ### Compute AF-$\beta_{LM}$s

print(f'\n* Calculating AF-BLMs...')

# data.afblmMatEfit(data = None)  # OK
BetaNormX, basis = data.afblmMatEfit()  # OK, uses default polarizations & ADMs as set in data['subset']
# BetaNormX, basis = data.afblmMatEfit(ADM = data.data['subset']['ADM'])  # OK, but currently using default polarizations
# BetaNormX, basis = data.afblmMatEfit(ADM = data.data['subset']['ADM'], pol = data.data['pol']['pol'].sel(Labels=['x']))
# BetaNormX, basis = data.afblmMatEfit(ADM = data.data['subset']['ADM'], pol = data.data['pol']['pol'].sel(Labels=['x','y']))  # This fails for a single label...?
# BetaNormX, basis = data.afblmMatEfit(RX=data.data['pol']['pol'])  # This currently fails, need to check for consistency in ep.sphCalc.WDcalc()
                                                                    # - looks like set values and inputs are not consistent in this case? Not passing angs correctly, or overriding?
                                                                    # - See also recently-added sfError flag, which may cause additional problems.


# ### AF-$\beta_{LM}$s

# The returned objects contain the $\beta_{LM}$ parameters as an Xarray...

# BetaNormX



# ## Fitting the data
#
# In order to fit data, and extract matrix elements from an experimental case, we'll use the [lmfit library](https://lmfit.github.io/lmfit-py/intro.html). This wraps core Scipy fitting routines with additional objects and methods, and is further wrapped for this specific class of problems in `pemtkFit` class we're using here.

# ### Set the data to fit
#
# Here we'll use the values calculated above as our test data. This currently needs to be set as `self.data['subset']['AFBLM']` for fitting.

# +
# data.data['subset']['AFBLM'] = BetaNormX  # Set manually

data.setData('sim', BetaNormX)  # Set simulated data to master structure as "sim"
data.setSubset('sim','AFBLM')   # Set to 'subset' to use for fitting.

# -

# Set basis functions
data.basis = basis

# ### Setting up the fit parameters
#
# In this case, we can work from the existing matrix elements to speed up parameter creation, although in practice this may need to be approached ab initio - nonetheless, the method will be the same, and the ab initio case detailed later.

# Input set, as defined earlier
# data.data['subset']['matE'].pd

print('\n*Setting  up fit parameters (with constraints)...')
# data.setMatEFit()  # Need to fix self.subset usage
# data.setMatEFit(data.data['subset']['matE'])  #, Eke=1.1) # Some hard-coded things to fix here! Now roughly working.

# With constraints
# Set param constraints as dict
# paramsCons = {}
# paramsCons['m_PU_SG_PU_1_n1_1'] = 'm_PU_SG_PU_1_1_n1'
# paramsCons['p_PU_SG_PU_1_n1_1'] = 'p_PU_SG_PU_1_1_n1'
#
# paramsCons['m_PU_SG_PU_3_n1_1'] = 'm_PU_SG_PU_3_1_n1'
# paramsCons['p_PU_SG_PU_3_n1_1'] = 'p_PU_SG_PU_3_1_n1'
#
# data.setMatEFit(paramsCons = paramsCons)

# With auto setting (from existing matrix elements)
data.setMatEFit()

print('\n\n*** Setup demo fitting workspace OK.')

# This sets `self.params` from the matrix elements, which are a set of (real) parameters for lmfit, as [a Parameters object](https://lmfit.github.io/lmfit-py/parameters.html).
#
# Note that:
#
# - The input matrix elements are converted to magnitude-phase form, hence there are twice the number as the input array, and labelled `m` or `p` accordingly, along with a name based on the full set of QNs/indexes set.
# - One phase is set to `vary=False`, which defines a reference phase. This defaults to the first phase item.
# - Min and max values are defined, by default the ranges are 1e-4<mag<5, -pi<phase<pi.
# - No relationships between the parameters are set by default (apart from the single fixed phase), but can be set manually, [see section below](http://127.0.0.1:8888/lab/workspaces/pemtk#Setting-parameter-relations).
