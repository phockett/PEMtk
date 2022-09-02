---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.2.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Basic PEMtk fitting class demo
01/09/22


Outline of this notebook:

- Load required packages.
- Setup pemtkFit object.
   - Set various parameters, either from existing data or new values.
   - Data is handled as a set of dictionaries within the class, `self.data[key][dataType]`, where `key` is an arbitrary label for, e.g. a specific experiment, calculation etc, and `dataType` contains a set of values, parameters etc. (Should become clear below!)
   - Methods operate on all `self.data` items in general, with some special cases: `self.data['subset']` contains data to be used in fitting.
- Simulate data
   - Use ePSproc to simulated aligned-frame measurements.
- Fit data (serial version)

For batch and parallel fitting see [the multi-fit demo notebook](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_demo_multi-fit_tests_130621-para.html)


### Prerequisities

- Working installation of [ePSproc](https://github.com/phockett/ePSproc) + [PEMtk](https://github.com/phockett/PEMtk) (or local copies of the Git repos, which can be pointed at for setup below).
- Test/demo data, from [ePSproc Github repo](https://github.com/phockett/ePSproc/tree/master/data/photoionization/).


## Setup


### Imports


A few standard imports...

```python
import sys
import os
from pathlib import Path
# import numpy as np
# import epsproc as ep
# import xarray as xr

from datetime import datetime as dt
timeString = dt.now()
```

And local module imports. This should work either for installed versions (e.g. via `pip install`), or for test code via setting the base path below to point at your local copies.

```python
# For module testing, include path to module here, otherwise use global installation
if sys.platform == "win32":
    modPath = Path(r'D:\code\github')  # Win test machine
    winFlag = True
else:
    modPath = Path(r'/home/femtolab/github')  # Linux test machine
    winFlag = False
    
# Append to sys path
sys.path.append((modPath/'ePSproc').as_posix())
sys.path.append((modPath/'PEMtk').as_posix())
```

```python
# ePSproc
import epsproc as ep

# Set data path
# Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
epDemoDataPath = Path(ep.__path__[0]).parent/'data'
```

```python
# PEMtk
# import pemtk as pm
# from pemtk.data.dataClasses import dataClass

# Import fitting class
from pemtk.fit.fitClass import pemtkFit
```

```python
# Set HTML output style for Xarray in notebooks (optional), may also depend on version of Jupyter notebook or lab, or Xr
# See http://xarray.pydata.org/en/stable/generated/xarray.set_options.html
# if isnotebook():
# xr.set_options(display_style = 'html')
```

```python
# Set some plot options
ep.plot.hvPlotters.setPlotters()
```

### Set & load parameters

There are a few things that need to be configured for a given case...

- Matrix elements (or (l,m,mu) indicies) to use.
- Alignment distribution (ADMs).
- Polarization geometry.

In this demo, a real case will first be simulated with computational values, and then used to test the fitting routines.

(TODO: demo from scratch without known matrix elements.)


#### Matrix elements

For fit testing, start with computational values for the matrix elements. These will be used to simulate data, and also to provide a list of parameters to fit later.

```python
# Set for ePSproc test data, available from https://github.com/phockett/ePSproc/tree/master/data
# Note this is set here from ep.__path__, but may not be correct in all cases.

# Multiorb data
dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
```

```python
data = pemtkFit(fileBase = dataPath, verbose = 1)
```

```python
# Read data files
data.scanFiles()
data.jobsSummary()
```

#### Alignment distribution moments (ADMs)

The class [wraps ep.setADMs()](https://epsproc.readthedocs.io/en/dev/modules/epsproc.sphCalc.html#epsproc.sphCalc.setADMs). This returns an isotropic distribution by default, or values can be set explicitly from a list. Values are set in `self.data['ADM']`.

Note: if this is not set, the default value will be used, which is likely not very useful for the fit!

```python
# Default case
data.setADMs()
# data.ADM['ADMX']
data.data['ADM']['ADM']
```

```python
# Load time-dependent ADMs for N2 case
# Adapted from ePSproc_AFBLM_testing_010519_300719.m

from scipy.io import loadmat
ADMdataFile = os.path.join(epDemoDataPath, 'alignment', 'N2_ADM_VM_290816.mat')
ADMs = loadmat(ADMdataFile)

# Set tOffset for calcs, 3.76ps!!!
# This is because this is 2-pulse case, and will set t=0 to 2nd pulse (and matches defn. in N2 experimental paper)
tOffset = -3.76
ADMs['time'] = ADMs['time'] + tOffset

data.setADMs(ADMs = ADMs['ADM'], t=ADMs['time'].squeeze(), KQSLabels = ADMs['ADMlist'], addS = True)
data.data['ADM']['ADM']
```

```python
# The ADMplot routine will show a basic line plot, note it needs keys = 'ADM' in the current implementation (otherwise will loop over all keys)
data.ADMplot(keys = 'ADM')
```

```python
# lmPlot is also handy, it minimally needs the correct keys, dataType and xDim specified

# data.lmPlot(keys = 'ADM', dataType = 'ADM', xDim = 't')  # Minimal call

data.lmPlot(keys = 'ADM', dataType = 'ADM', xDim = 't', fillna = True, logFlag = False, cmap = 'vlag')  # Set some additional options
```

### Polarisation geometry/ies

This wraps [ep.setPolGeoms](https://epsproc.readthedocs.io/en/dev/modules/epsproc.sphCalc.html#epsproc.sphCalc.setPolGeoms). This defaults to (x,y,z) polarization geometries. Values are set in `self.data['pol']`.

Note: if this is not set, the default value will be used, which is likely not very useful for the fit!


```python
data.setPolGeoms()
data.data['pol']['pol']
```

```python
# # data.setPolGeoms(eulerAngs = [[0,0,0]], labels = ['z'])
# data.setPolGeoms(eulerAngs = [0,0,0], labels = 'z')
# data.data['pol']['pol']  #.swap_dims({'Euler':'Labels'})
```

```python
# data.data['pol']['pol'] = data.data['pol']['pol'].swap_dims({'Euler':'Labels'})
# data.selOpts['pol'] = {'inds': {'Labels': 'z'}}
# data.setSubset(dataKey = 'pol', dataType = 'pol')
```

### Subselect data

Currently handled in the class by setting `self.selOpts`, this allows for simple reuse of settings as required. Subselected data is set to `self.data['subset'][dataType]`, and is the data the fitting routine will use.

```python
# Settings for type subselection are in selOpts[dataType]

# E.g. Matrix element sub-selection
data.selOpts['matE'] = {'thres': 0.01, 'inds': {'Type':'L', 'Eke':1.1}}
data.setSubset(dataKey = 'orb5', dataType = 'matE')  # Subselect from 'orb5' dataset, matrix elements

# Show subselected data
data.data['subset']['matE']
```

```python
# Tabulate the matrix elements
# Not showing as nice table for singleton case - pd.series vs. dataframe?
data.matEtoPD(keys = 'subset', xDim = 'Sym', drop=False)

# data.data['subset']['matE'].attrs['pd'] # PD set here
```

```python
# And for the polarisation geometries...
data.selOpts['pol'] = {'inds': {'Labels': 'z'}}
data.setSubset(dataKey = 'pol', dataType = 'pol')
```

```python
# And for the ADMs...

data.selOpts['ADM'] = {}   #{'thres': 0.01, 'inds': {'Type':'L', 'Eke':1.1}}
data.setSubset(dataKey = 'ADM', dataType = 'ADM', sliceParams = {'t':[4, 5, 4]}) 
data.ADMplot(keys = 'subset')
```

```python
# Cusomise plot with return...
# NOT YET IMPLEMENTED
# pltObj = data.ADMplot(keys = 'subset')
# pltObj
```

```python
# Plot from Xarray vs. full dataset
# data.data['subset']['ADM'].where(ADMX['K']>0).real.squeeze().plot.line(x='t');
data.data['subset']['ADM'].real.squeeze().plot.line(x='t', marker = 'x', linestyle='dashed');
data.data['ADM']['ADM'].real.squeeze().plot.line(x='t');
```

## Compute AF-$\beta_{LM}$ and simulate data

With all the components set, some observables can be calculated. For testing, we'll also use this to simulate an experiemental trace...

Here we'll use `self.afblmMatEfit()`, which is also the main fitting routine, and essentially wraps `epsproc.afblmXprod()` to compute AF-$\beta_{LM}$s (for more details, see the [ePSproc method development docs](https://epsproc.readthedocs.io/en/dev/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html)).

If called without reference data, the method returns computed AF-$\beta_{LM}$s based on the input subsets already created, and also a set of (product) basis functions generated - these can be examined to get a feel for the sensitivity of the geometric part of the problem, and will also be used in fitting to limit repetitive computation.


### Compute AF-$\beta_{LM}$s

```python
# data.afblmMatEfit(data = None)  # OK
BetaNormX, basis = data.afblmMatEfit()  # OK, uses default polarizations & ADMs as set in data['subset']
# BetaNormX, basis = data.afblmMatEfit(ADM = data.data['subset']['ADM'])  # OK, but currently using default polarizations
# BetaNormX, basis = data.afblmMatEfit(ADM = data.data['subset']['ADM'], pol = data.data['pol']['pol'].sel(Labels=['x']))
# BetaNormX, basis = data.afblmMatEfit(ADM = data.data['subset']['ADM'], pol = data.data['pol']['pol'].sel(Labels=['x','y']))  # This fails for a single label...?
# BetaNormX, basis = data.afblmMatEfit(RX=data.data['pol']['pol'])  # This currently fails, need to check for consistency in ep.sphCalc.WDcalc()
                                                                    # - looks like set values and inputs are not consistent in this case? Not passing angs correctly, or overriding?
                                                                    # - See also recently-added sfError flag, which may cause additional problems.

```

### AF-$\beta_{LM}$s


The returned objects contain the $\beta_{LM}$ parameters as an Xarray...

```python
BetaNormX
```

```python
# ep.BLMplot(BetaNormX, xDim = 't')  # SIGH, not working - issue with Euler/Labels?
# BetaNormX.sel(Labels='z').real.squeeze().plot.line(x='t');
# ep.lmPlot(BetaNormX.sel(Labels='z'), xDim='t', SFflag=False);  #, cmap='vlag');
ep.lmPlot(BetaNormX, xDim='t', SFflag=False);
# data.lmPlot()
```

```python
# Line-plot with Xarray/Matplotlib
# Note there is no filtering here, so this includes some invalid and null terms
BetaNormX.sel(Labels='A').real.squeeze().plot.line(x='t');
```

... and the basis sets as a dictionary.

```python
basis.keys()
```

Note that the basis sets here will may be useful for deeper insight into the physics, and fitting routines, and are explored in a separate notebook.


## Fitting the data

In order to fit data, and extract matrix elements from an experimental case, we'll use the [lmfit library](https://lmfit.github.io/lmfit-py/intro.html). This wraps core Scipy fitting routines with additional objects and methods, and is further wrapped for this specific class of problems in `pemtkFit` class we're using here.


### Set the data to fit

Here we'll use the values calculated above as our test data. This currently needs to be set as `self.data['subset']['AFBLM']` for fitting.

```python
# data.data['subset']['AFBLM'] = BetaNormX  # Set manually

data.setData('sim', BetaNormX)  # Set simulated data to master structure as "sim"
data.setSubset('sim','AFBLM')   # Set to 'subset' to use for fitting.

```

```python
# Set basis functions
data.basis = basis
```

### Setting up the fit parameters

In this case, we can work from the existing matrix elements to speed up parameter creation, although in practice this may need to be approached ab initio - nonetheless, the method will be the same, and the ab initio case detailed later.

```python
# Input set, as defined earlier
data.data['subset']['matE'].pd
```

```python
# data.setMatEFit()  # Need to fix self.subset usage
data.setMatEFit(data.data['subset']['matE'])  #, Eke=1.1) # Some hard-coded things to fix here! Now roughly working.
```

This sets `self.params` from the matrix elements, which are a set of (real) parameters for lmfit, as [a Parameters object](https://lmfit.github.io/lmfit-py/parameters.html). 

Note that: 

- The input matrix elements are converted to magnitude-phase form, hence there are twice the number as the input array, and labelled `m` or `p` accordingly, along with a name based on the full set of QNs/indexes set.
- One phase is set to `vary=False`, which defines a reference phase. This defaults to the first phase item.
- Min and max values are defined, by default the ranges are 1e-4<mag<5, -pi<phase<pi.
- Relationships between the parameters are set by default, but can be set manually, [see section below](#Setting-parameter-relations), or pass `paramsCons=None` to skip.

```python
data.params
```

### Running a fit...

With the parameters and data set, just call `self.fit()`!

Statistics and outputs are handled by lmfit, which includes uncertainty estimates and correlations in the fitted parameters.

```python
data.fit()
```

```python
# Check fit outputs - self.result shows results from the last fit
data.result
```

Results vs. data can be (crudely) plotted with `self.BLMfitPlot()` (better plotting routines to follow!).

This will plot all results by default, vs. the `subset` data used for the fitting routine inputs.

```python
# Plot data subset (--x) plus fit (solid lines)
data.BLMfitPlot()
```

```python
# Fit results are currently added to the main data dict by an index number
data.data.keys()
```

```python
# The current fit index is set in self.fitInd
data.fitInd
```

```python
# Full results are available from the main data structure
data.data[data.fitInd].keys()
```

```python
# Plot results with lmPlot wrapper - this also defaults to show data subset + most recent fit results
data.lmPlotFit()
```

### Fitting with randomised parameter inputs

Here we might expect some variation in the results, depending on various properites (ionizing channel, dataset size etc.), and also for the fitting to take a bit longer (see benchmarks later).

TODO: 

- More careful analysis routines here, run vs. number of input points & test fiedelity.
- Parallelize.
- Variation over runs/general statistical analysis.

```python
# Basic randomize routine, [0,1] interval
data.randomizeParams()
```

```python
data.fit()
# data.data[data.fitInd]['results']
data.result
```

```python
# Note that if keys is not set, BLMfitPlot will show all fit run results.
data.BLMfitPlot()
```

### Setting parameter relations/constraints

If we know that some of the matrix elements (parameters) are related, this can be set [using contraints on the parameters](https://lmfit.github.io/lmfit-py/constraints.html).

In this test case, we know some terms are equal... this should speed up the fitting, and also improve fiedelity.

UPDATE Aug 2022: note that constraints are now set automatically by `self.setMatEFit()` if working from known matrix elements. This can be overridden and/or set manually if desired.

```python
# With auto constraints - this is the default case
data.setMatEFit(paramsCons = 'auto')
```

```python
# With NO constraints set (except a reference phase)
data.setMatEFit(paramsCons = None)
```

```python
# With manual constraints
# Set param constraints as dict
paramsCons = {}
paramsCons['m_PU_SG_PU_1_n1_1_1'] = 'm_PU_SG_PU_1_1_n1_1'
paramsCons['p_PU_SG_PU_1_n1_1_1'] = 'p_PU_SG_PU_1_1_n1_1'

paramsCons['m_PU_SG_PU_3_n1_1_1'] = 'm_PU_SG_PU_3_1_n1_1'
paramsCons['p_PU_SG_PU_3_n1_1_1'] = 'p_PU_SG_PU_3_1_n1_1'

# Missing settings will generate an error message
paramsCons['test'] = 'p_PU_SG_PU_3_1_n1_1'

data.setMatEFit(paramsCons = paramsCons)
```

```python
data.randomizeParams()
data.fit()
data.result
```

```python
# Note that if keys is not set, BLMfitPlot will show only most recent fit run results.
data.BLMfitPlot()
```

### Quick benchmarks

<!-- #region -->
#### Quick benchmark - basic fit (no constraints)

Timing for the test fit (next cell) - note this is currently just running with defaults, single core only, and we're not checking for good fiedelity here. Fits take between 1 and 12 mins, approximately.

Results may vary depending on the inputs... for the current test case:

- 1st test
    The slowest run took 12.88 times longer than the fastest. This could mean that an intermediate result is being cached.
    50.7 s ± 1min 2s per loop (mean ± std. dev. of 7 runs, 1 loop each)
- 2nd test
    The slowest run took 105.38 times longer than the fastest. This could mean that an intermediate result is being cached.
    4min 13s ± 7min 6s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    
    
TODO:

- More careful testing & benchmarks.
- Timing vs. input dataset size.
- Fitting statistics over fits (fidelity vs. fit vs. dataset size etc.).
<!-- #endregion -->

```python
# %%timeit

# # With NO constraints set (except a reference phase)
# data.setMatEFit(paramsCons = None)

# data.randomizeParams()
# data.fit()
```

#### Quick benchmark - constrained fit

Timing for the test fit (next cell, same [method as earlier](http://127.0.0.1:8888/lab/workspaces/pemtk#Quick-benchmark---basic-fit)) - note this is currently just running with defaults, single core only, and we're not checking for good fiedelity here. Note, also, that these results are in the ~10s range, compared to ~many minutes for the unconstrained case.

- 1st test
    10.7 s ± 2.86 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
- 2nd test 7.9 s ± 3.5 s per loop (mean ± std. dev. of 7 runs, 1 loop each)


```python
%%timeit

data.randomizeParams()
data.fit()
```

```python
# Plot multiple fit sets
# Note that if keys = 'all' is set, BLMfitPlot will show ALL fit run results.
# TODO: use Seaborn/HV for better plotting here!
data.BLMfitPlot(keys = 'all')
```

## Versions

```python
import scooby
scooby.Report(additional=['epsproc', 'pemtk', 'xarray', 'jupyter'])
```

```python
# Check current Git commit for local ePSproc version
!git -C {Path(ep.__file__).parent} branch
!git -C {Path(ep.__file__).parent} log --format="%H" -n 1
```

```python
# Check current remote commits
!git ls-remote --heads https://github.com/phockett/ePSproc
# !git ls-remote --heads git://github.com/phockett/epsman
```

```python
# Check current Git commit for local PEMtk version
import pemtk
!git -C {Path(pemtk.__file__).parent} branch
!git -C {Path(pemtk.__file__).parent} log --format="%H" -n 1
```

```python
# Check current remote commits
!git ls-remote --heads https://github.com/phockett/PEMtk
# !git ls-remote --heads git://github.com/phockett/epsman
```

```python

```
