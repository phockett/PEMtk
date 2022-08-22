# PEMtk wrappers for parallel processing for fitting
#
# 13/09/21  v1  Basic wrapper with XYZPY implemented,
#
# Initial dev code: see http://127.0.0.1:8888/lab/tree/dev/PEMtk/fitting/multiFit_tests_and_parallel/PEMtk_fitting_dev_multiproc_070921.ipynb
#

# Currently using xyzpy for parallel stuff.
try:
    import xyzpy as xyz
except:
    print("***xyzpy not found, parallel functions not available.")

import multiprocessing as mp

def multiFit(self, nRange = [0,10],
            parallel = True, num_workers = None,
            randomizeParams = True, seedParams = None):
    """
    Basic wrapper for pemtk.fitClass.fit() for multiprocess execution.

    Run a batch of fits in parallel, and return results to main class structure.


    Parameters
    ----------
    nRange : list
        Fit indexers.
        Set [nStart, nStop], full run will be set as list(range(nRange[0],nRange[1])).
        TODO: more flexibility here, and auto.

    parallel : bool, default = True
        Run fit jobs in parallel?

    num_workers : int, default = None
        Number of cores to use if parallel job.
        Currently set to default to ~90% of mp.cpu_count()

    randomizeParams : bool, default = True
        Randomize seed parameters per fit?

    seedParams : int, default = None
        NOT IMPLEMENTED, but will provide an option to seed fits from a previous result.



    """

    #*** Settings
    if num_workers is None:
        num_workers =  round(mp.cpu_count() * 0.9)   # Default to ~90% cores

    if self.verbose['main']:
        print("Number of processors: ", mp.cpu_count(), "\nRunning pool on: ", num_workers)


    #*** Setup XYZPY runner
    combos = {'n': list(range(nRange[0],nRange[1]))}
    data = {'data':self}   # Pass data - may be better to set fitPara as class method? Not sure if that will work with wrapper.

    # With runner class
    r = xyz.Runner(fitPara, var_names=['results'], resources = data, constants = {'randomizeParams':randomizeParams, 'seedParams':seedParams})  # Use constants to pass additional function args.

    #*** Run fits
    outputs = r.run_combos(combos, parallel = parallel, num_workers = num_workers)  # This returns an Xarray Dataset.

    #*** Sort results back to main data structure
    # Restack from output > dict form
    # Bit hacky, but works OK.
    for n in outputs.n:
        self.data[n.item()] = outputs.results.sel({'n':n}).item()

    # Update self.fitInd 
    fInd, fitInds = self._getFitInds()
    self.fitInd = fInd


def fitPara(data = None, n = None, randomizeParams = True, seedParams = None):
    """
    Wrap self.fit() for XYZPY runner.
    """

    if randomizeParams:
        data.randomizeParams()  # This may not work as is, since it might accidentally overwrite in-use params?
                                # Comment out to give fast test run only
    if seedParams is not None:
        # data.params = data.data[seedParams]['results']    # Should be something like this, just need to pull correct vars.
        print('seedParams not yet implemented.')


    data.fit(fitInd = n)    # This may also need some attention if pointers are shared.

    return data.data[n]  # OK
