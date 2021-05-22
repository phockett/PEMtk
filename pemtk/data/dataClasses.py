"""
PEMtk base data classes

Main data structures & types for experimental & computational data, latter based on ePSproc data types.

NOTE - this requires ePSproc, local import currently set in top level __init__.py

13/10/20    v1

"""

# Multijob class dev code
from epsproc.classes.multiJob import ePSmultiJob
from epsproc.classes.base import ePSbase

# For basic testing, just wrap ePSproc multiJob class.
# NOTE - this requires ePSproc, local import currently set in top level __init__.py
class dataClass(ePSbase):
    """
    PEMtk data class stub - currently just a thin wrapper to ePSproc.classes.ePSbase for dev/testing.

    """

    pass
