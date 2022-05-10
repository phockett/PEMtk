# PEMtk fitting stats routines
#
# Helper functions for sampling & stats.
#
# 04/05/22  v1  Adding some basic methods for sampling & bootstrapping.
#
# TODO:
#   - Sampling methods for Xarrays: https://github.com/smartass101/xr-random
#   - Bootstrapping with Arch library: https://arch.readthedocs.io/en/latest/bootstrap/bootstrap.html
#

import numpy as np

def setPoissWeights(lam, wShape):
    """
    Poissionian weights.

    Set using https://numpy.org/devdocs/reference/random/generated/numpy.random.Generator.poisson.html

    Pass (lambda, shape)

    If int or float, set weights = rng.poisson(weights, data.shape).
    Note this operates per data-point, not per dimension.

    - If Xarray or np.array, use directly - must match size of data along key dimension, e.g. passing weights = rng.poisson(weights, data.t.size) will generate a distribution along the t-dimension.

    """

    rng = np.random.default_rng()
    weights = rng.poisson(lam, wShape)

    return weights
