"""
Utility functions for symmetrized harmonics routines.

23/03/22

"""

import numpy as np

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
