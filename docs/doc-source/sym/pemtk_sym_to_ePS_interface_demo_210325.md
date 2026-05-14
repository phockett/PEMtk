---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Symmetrized harmonics demo pt IV: photoionization selection rules for ePolyScat inputs
21/03/25

This notebook demonstrates functionality for working with symmetry, making use of [libmsym](https://github.com/mcodev31/libmsym) for the core symmetry functionality.

The focus here is on symmetry allowed ionization channels, and setting appropriate variables for [ePolyScat (ePS) input files](https://epolyscat.droppages.com/).

For more background, see the [Symmetriezed Haromnics demo](https://pemtk.readthedocs.io/en/latest/sym/pemtk_symHarm_demo_160322_tidy.html).

For more on point groups and symmetry, see http://www.gernot-katzers-spice-pages.com/character_tables/index.html, which includes full character tables and associated data.

+++

## Generate basis

```{code-cell} ipython3
# Import class and generate basis
from pemtk.sym.symHarm import symHarm

sym = 'D10h'
symBasis = symHarm(PG = sym, lmax=4)
```

## Set ionization channel & compute symmetries

For more options and discussion, see [the demo pt III](https://pemtk.readthedocs.io/en/latest/sym/working_with_symmetry_PEMtk_notes_180423-v210423_tidy.html#Mapping-allowed-terms-to-photoionization-matrix-elements).

```{code-cell} ipython3
# Specify full terms - this will indicate if a given selection fulfils the dipole selection rules for the dipole terms in the given point group.
sNeutral = 'A1g'
sIon = 'E1u'
sContinuum = 'A1g'

# The `directProductDipole` method outputs results with the dipole symmetries included
symBasis.directProductDipole([sNeutral, sIon, sContinuum])
symBasis.dipole['dipoleProducts']
```

For this case we see that the set of states `[sNeutral, sIon, sContinuum]` is allowed for the $E_{1g}(x,y)$ dipole component.

```{code-cell} ipython3
# For cases where the continuum symmetry is not known/specified, the `directProductContinuum` method can be used, and will output results looped over all irreps
symBasis.directProductContinuum([sNeutral, sIon])
```

```{code-cell} ipython3
# For ePS, we also need to define the scattering symmetry, defined as scatSym = continuum x ion
# Run self.scatSym to compute this
symBasis.scatSym(sIon)
```

## Output to ePS/ePSman format

Once computed, the necessary terms are output to a list of symmety tuples `(Total, Cont)`, where "Total" is the scattering symmetry, and "Cont" the corresponding continuum (photoelectron) symmetry.

```{code-cell} ipython3
symBasis.continuum['allowed']['ePSSymList']
```

Note some other lower-level methods are available, in particular `self.toePSman()` and `self.toePSmanPD()` for conversions from list or Pandas dataframes.

TODO: explicitly convert list to ePS symbols (currently this can be done in ePSman).

+++

## Additional notes

+++

### Multiple unpaired electrons

In some cases the initial or final state may require multiple characters, this can be handled by setting the states as lists and computing the relevant direct product(s).

```{code-cell} ipython3
# Setting the symmetry for the neutral and ion allows direct products to be computed, 
# and allowed terms to be determined.

sNeutral = ['B2u','B1g']
sIon = 'B2u'
```

```{code-cell} ipython3
# Check total neutral state character
symBasis.directProduct(sNeutral)
```

```{code-cell} ipython3
# Compute continuum results for multiple neutral UPEs - simply pass as list
symBasis.directProductContinuum([*sNeutral, sIon])

# Results are pushed to self.continuum, in dictionary and Pandas DataFrame formats, 
# and can be manipulated using standard functionality.
# The subset of allowed values are also set to a separate DataFrame and list.
symBasis.continuum['allowed']['PD']
```

## Versions

```{code-cell} ipython3
import scooby
scooby.Report(additional=['epsproc', 'pemtk', 'holoviews', 'hvplot', 'xarray', 'matplotlib', 'bokeh', 'pyshtools'])
```

```{code-cell} ipython3
# Check current Git commit for local ePSproc version
import pemtk
from pathlib import Path
!git -C {Path(pemtk.__file__).parent} branch
!git -C {Path(pemtk.__file__).parent} log --format="%H" -n 1
```

```{code-cell} ipython3
# Check current remote commits
!git ls-remote --heads https://github.com/phockett/pemtk
```
