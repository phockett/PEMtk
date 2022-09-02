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

# PEMtk fitting setup & batch run demo
01/09/22

Outline of this notebook:

- Use `setup_fit_demo.py` script to load data and setup fitting environment.
- Run a batch of fits...
    
    - (a) [run batch](#(a)-Run-a-batch)
        - Serial with `self.fit()` method.
        - [Parallel fitting routine](#(2)-parallel-execution) with `self.multiFit()` method. This uses [XYZpy library](https://xyzpy.readthedocs.io/en/latest/) for quick parallelization.
        
    - (b) load a batch
    
- Exploring the results... see [the "analysis" notebook](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy.html) for more.

---

Versions

- 06/06/21  v1
- 22/08/22  v2
- 01/09/22  v2b

v2b updated version including:

- Data IO with `self.writeFitData()` method.
- fitting now supports additional backends, these are set in `self.fitOpts`. For a list of available backends run `self.backends()`. Currently this implements support for AF or MF fitting routines.


## Setup

```python
# Import & set paths
import pemtk
from pemtk.fit.fitClass import pemtkFit

from pathlib import Path

# Path for demo script
demoPath = Path(pemtk.__file__).parent.parent/Path('demos','fitting')

# Run demo script to configure workspace
%run {demoPath/"setup_fit_demo.py"}
```

## Run fits

For this notebook, we'll either (a) run a batch of fits or (b) load sample data.

With the current codebase, running multiple fits will default to using the same basis set, and output results sequentially to the main `self.data` dictionary.

Update 20/08/22: fitting now supports additional backends, these are set in `self.fitOpts`. For a list of available backends run `self.backends()`. Currently this implements support for AF or MF fitting routines.

Benchmarks for an [AMD Threadripper 2950X (16 core) system](#Versions).


### (a) Run a batch

#### (1) serial execution

Either:

- Manually with a loop.
- With `self.multiFit()` method, although this is optimised for parallel execution (see below).

```python
import time

start = time.time()

# Maual execution
for n in range(0,100):
    data.randomizeParams()
    data.fit()
    
end = time.time()
print((end - start)/60)
    
# Or run with self.multiFit(parallel = False)
# data.multiFit(nRange = [0,100], parallel = False)
```

```python
# We now have 100 fit results
data.data.keys()
```

```python
# Quick data dump
# TODO: better save routine (json/h5).
# Now wrapped in self.writeFitData(), see https://github.com/phockett/PEMtk/issues/6 and docs at https://epsproc.readthedocs.io/en/dev/dataStructures/ePSproc_dataStructures_IO_demo_280622.html

# import pickle
# with open('dataDump_100fitTests_130621.pickle', 'wb') as handle:
#     pickle.dump(data.data, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Data IO with self.writeFitData()
# This will default to the working dir and set a data-stamped file name if nothing is passed.
data.writeFitData()
```

#### (2) parallel execution

Updated version including parallel fitting routine with `self.multiFit()` method. 

This currently uses the [XYZpy library](https://xyzpy.readthedocs.io/en/latest/) for quick parallelization, although there is some additional setup overhead in the currently implementation due to class init per fit batch. The default aims to set ~90% CPU usage, based on core-count.

```python
data.multiFit(nRange = [0,100])
```

### (b) Load a batch of fit runs

Load sample data for analysis instead of running fits. Note this can be run minimally *without* the full setup routines above, using the commented-out cell below to init a blank object.

(The [demo file(s) are available in demos/fitting](https://github.com/phockett/PEMtk/tree/master/demos/fitting).)

```python
# If running from scratch, create a blank object first
# # Init blank object
from pemtk.fit.fitClass import pemtkFit
data = pemtkFit()
```

```python
# Load sample dataset
# Full path to the file may be required here, in demos/fitting
# import pickle

# with open('dataDump_100fitTests_10t_randPhase_130621.pickle', 'rb') as handle:
#     data.data = pickle.load(handle)
    
# Data IO with self.writeFitData()
# This will default to the working dir
data.loadFitData('dataDump_100fitTests_10t_randPhase_130621.pickle')
```

```python
data.data.keys()
```

## Exploring a fit result

Each result contains a set of fit results:

```python
nFit = 11
data.data[nFit].keys()
```

Here 'results' is an [lmFit object](https://lmfit.github.io/lmfit-py/intro.html), containing various outputs, including the final paramter set and fit statistics, which can be inspected directly. (See the [basic demo notebook for more.](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html))

```python
data.data[nFit]['results']
```

The best fit results are set in an Xarray, keyed by `AFBLM`.

```python
data.data[nFit]['AFBLM']
```

For further analysis & batch results, see [the "analysis" notebook](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy.html).


## Versions

```python
import scooby
scooby.Report(additional=['epsproc', 'pemtk', 'xarray', 'jupyter'])
```

```python
# Check current Git commit for local ePSproc version
from pathlib import Path
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
from pathlib import Path
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
