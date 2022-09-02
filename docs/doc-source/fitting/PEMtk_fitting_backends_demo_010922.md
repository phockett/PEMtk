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

# Fitting model functions (backends)
01/09/22

As of August 2022, both AF and MF model functions are supported for fitting, and additional models may now readily be added via the `self.backends()` routine. This notebook briefly demos the functionality.

(For a general fitting demo, see the [basic demo notebook](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html).)


## Setup

```python
# Import fit class
from pemtk.fit.fitClass import pemtkFit

# Init new fit class object
data = pemtkFit(verbose = 2)
```

## Default model backend


The default backend is set at class init, and the name set in `self.backend`

```python
data.backend
```

For use in fitting, the function corresponding to this name is set in `self.fitOpts['backend']`.

```python
data.fitOpts['backend']
```

## Available models

Call `self.backends()` to get a dictionary of available functions and handles.

```python
data.backends()
```

To change the current model, call .backends() with a supported function name from the list above to set.

```python
data.backends('mfblmXprod');
```

```python
data.backend
```

```python
data.fitOpts['backend']
```

```python
data.backends('afblmXprod');
```

```python
data.backend
```

## Custom backends

To use alternative backends, just set the function handle to `self.fitOpts['backend']`, this can also be set at class creation.

```python
# Init new fit class object
import numpy as np
data2 = pemtkFit(verbose = 2, backend=np.abs)
```

```python
data2.backend
```

```python
data2.fitOpts['backend']
```

Note, however, that this function is currently used by the `self.afblmMatEfit`, and must return an Xarray of BLM parameters, and a set of basis functions - see [the source code for more details](https://pemtk.readthedocs.io/en/latest/_modules/pemtk/fit/fitClass.html#pemtkFit.afblmMatEfit).

In future, more general wrappers may be added here.


## Versions

```python
import scooby
scooby.Report(additional=['epsproc', 'pemtk', 'xarray', 'jupyter'])
```

```python
# Check current Git commit for local ePSproc version
import epsproc as ep
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
import epsproc as ep
ep.__file__
```

```python

```
