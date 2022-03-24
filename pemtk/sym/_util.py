"""
Utility functions for symmetrized harmonics routines.

23/03/22

"""

import numpy as np

try:
    from epsproc.util.listFuncs import dataTypesList
except ImportError:
    print("***Missing epsproc, can't import data types.")


def prefixAxisTitles(p, prefix = '', keepTitle = True):
    """
    General get/set title for SHtools plot objects.

    Sets title to prefix, plus any existing text.
    Bit of a hack, but fixes issues with missing grid.plot(title) functionality in some versions.

    Required in testing for pyshtool v4.5, matplotlib v3.2.0.
    Not required for pyshtools v4.9.

    Parameters
    ----------

    p : grid.plot() returned object, usually (fig, [axes]), may be nested.

    prefix : str, optional, default = ''
        Prefix string to apply.
        If keepTitle = False this will be used as the title.

    keepTitle : bool, optional, default = True
        Keep existing title text?

    """

    # p usually (fig, [axes]), may be nested.
    for item in p:
#         print(item)

        # Axis object has set/get attribs
        if hasattr(item, 'get_title'):
            title = item.get_title()

            if title and prefix and keepTitle:
#                 print('Setting' + prefix + ', ' + title)
                item.set_title(prefix + ', ' + title)

            else:
                item.set_title(prefix)

#             print(f'Set {prefix + title}')
#             print(item)

        # Nested axis objects are np.ndarrays. Parse recursively.
        elif type(item) is np.ndarray:
#             print(type(item))
            prefixAxisTitles(item, prefix)  # Call recursively for nested axes (as ndarray)

def listPGs():
    """Return list of supported Point Groups (libmsym, https://github.com/mcodev31/libmsym)."""

    return ['Ci', 'Cs', 'Cnv', 'Dn', 'Dnh', 'Dnd', 'Td', 'O', 'Oh', 'I', 'Ih']



def toePSproc(coeffs, dimMap = {'C':'Cont','h':'it'}, dataType = 'BLM'):
    """
    Convert/conform Xarray of symmetrized harmonics to reference ePSproc data structures.

    This allows for quick use of general functions for matrix elements and BLM parameters (e.g. plotting), although may not work in all cases.

    Note this requires ep.util.listFuncs for definitions.

    Parameters
    ----------

    coeffs : dictionary
        Contains symmetrized harmonic coeffs, as per :py:func:`pemtk.symHarm.symHarm()` class definitions.

    dimMap : dictionary, default = {'C':'Cont','h':'it'}
        Any specific dim remapping required.
        The default case remaps to ePSproc labels 'Cont' and 'it'.

    dataType : str, default = 'BLM'
        Data type to conform to.
        See :py:func:`ep.dataTypesList()` for all supported types.
        For harmonics, 'BLM' or 'matE' are suggested.

    Returns
    -------
    XR : renamed & restacked DataArray


    """

    # dimMap = {'C':'Cont','h':'it'}  # Remap to ePS definitions - mu should also remap?
                                    # Not sure where to send it, should be h+mu for degen?
                                    # Or sum over mu?
    daTest = coeffs['XR'].copy()

    daTest = daTest.unstack()

    daTest = daTest.rename(dimMap)  # Remap existing dims

    # Set refDims from ep dataTypes function
    refDims = dataTypesList()[dataType]['def']  #(sType='unstacked')

    # Add dims as required
    for dim in refDims(sType='unstacked'):
        if not dim in daTest.dims:

            # Works OK, but get issues with Eke functionality later (should be int/float)
    #         daTest = daTest.expand_dims({dim: ['U']})

            # Get issues with lmplot for Nan or None case...
            # Maybe OK if dims are squeezed out first?
    #         daTest = daTest.expand_dims({dim: [np.nan]})
    #         daTest = daTest.expand_dims({dim: [None]})

            # Set for int vs. str labelled dims.
            # This works with existing lmplot routine OK
            if dim in ['Cont','Targ','Total','Type']:
                daTest = daTest.expand_dims({dim: ['U']})
            else:
                daTest = daTest.expand_dims({dim: [0]})
    #         daTest = daTest.expand_dims({dim: [0]})

    # Restack for output
    daTest = daTest.stack(refDims(sType='sDict'))

    daTest.attrs = coeffs['XR'].attrs  # Propagate attribs
    daTest.attrs['dataType'] = dataType

    return daTest
