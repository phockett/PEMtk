# IO functions for pemtk data
# For main matrix element routines etc. see ePSproc base class.
#
# 20/05/22  v1

# from datetime import datetime as dt
# import pickle

# Quick Pickle for all data in class.
# def pickleData(self, fName = None, outStem = ):
#
#     if fName is None:
#         timeString = dt.now()
#         fOut = f'{outStem}_n{n}_{timeString.strftime("%d%m%y_%H-%M-%S")}.pickle'
#
#     with open(fOut, 'wb') as handle:
#         pickle.dump(data.data, handle, protocol=pickle.HIGHEST_PROTOCOL)


import pickle
from pathlib import Path

from datetime import datetime as dt

import numpy as np

# dataPath = Path(r'/home/jovyan/work/pemtk_fitting_runs_April2022/ioTests')
# dataPath = Path(r'/home/jovyan/ioTestsLocal/data_dumps')  # FF Docker

# import os

def setTimeStampedFileName(outStem = None, n = None, ext = 'pickle', timeString = None):
    """
    Set unique filename as f'{outStem}_n{n}_{timeString.strftime("%d%m%y_%H-%M-%S")}.{ext}'

    Parameters
    ----------

    outStem : str, optional, default = None
        Stem for output file.
        If None, set to 'PEMtk_data_dump'

    n : int, optional, default = None
        Int index to include in file name.
        If None, this will be omitted.

    ext : str, optional, default = 'pickle'
        File ending.

    timeString : Datatime object, default = None
        Timestamp for the file.
        If None, current time will be used.

    """

    if outStem is None:
        outStem = 'PEMtk_data_dump'

    if timeString is None:
        timeString = dt.now()

    # Use n for multiple fit checkpointing - may want a more sophisticated routine for this however.
    if n is not None:
        fName = f'{outStem}_n{n}_{timeString.strftime("%d%m%y_%H-%M-%S")}.{ext}'
    else:
        fName = f'{outStem}_{timeString.strftime("%d%m%y_%H-%M-%S")}.{ext}'

    return fName


def writeFitData(self, dataPath = None, fName = None, outStem = None, n=None, fType = 'pickle', ext = None):
    """
    Dump fit data with various backends.



    """

    if dataPath is None:
        # dataPath = os.getcwd()  # OR
        dataPath = Path().absolute()

    if outStem is None:
        outStem = 'dataDump'

    if ext is None:
        ext = fType

    if fName is None:
        fName = setTimeStampedFileName(outStem = outStem, n = n, ext = ext)


    fOut = Path(dataPath,fName)

    if fType == 'pickle':
        # pickleData(self,fOut)  # NOTE - here only include self for testing!
        self._pickleData(fOut)

#     if self.verbose['main']:
#         print(f'Dumped data to {fOut} with pickle.')

    return fOut


# v1
# def pickleData(self, dataPath = None, fName = None, outStem = None, n=None):
#     """
#     Dump self.data to file using Pickle.

#     Usual Pickle caveats apply: not recommended for longevity, but handy for checkpoints and backup.

#     """

#     if dataPath is None:
#         # dataPath = os.getcwd()  # OR
#         dataPath = Path().absolute()

#     if outStem is None:
#         outStem = 'dataDump'

#     if fName is None:
#         fName = setTimeStampedFileName(outStem = outStem, n = n, ext = 'pickle')


#     fOut = Path(dataPath,fName)
#     with open(fOut, 'wb') as handle:
#         pickle.dump(data.data, handle, protocol=pickle.HIGHEST_PROTOCOL)

#     if self.verbose['main']:
#         print(f'Dumped data to {fOut} with pickle.')

# v2 - moved general path handling to wrapper
def _pickleData(self, fOut):
    """
    Dump self.data to file using Pickle.

    Usual Pickle caveats apply: not recommended for longevity, but handy for checkpoints and backup.

    """

    with open(fOut, 'wb') as handle:
        pickle.dump(self.data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if self.verbose['main']:
        print(f'Dumped self.data to {fOut} with pickle.')


def loadData(self, fList, dataPath = None, batch = None):
    """
    Load data dumps from a file or set of files (and stack).

    """

    if not isinstance(fList,list):
        fList = [fList]

    if dataPath is None:
        # dataPath = os.getcwd()  # OR
        dataPath = Path().absolute()

    fOffset = 0
    for ind,fileIn in enumerate(fList):
#         data.verbose['sub'] = 1

        if not Path(fileIn).exists():
            fileIn = Path(dataPath,fileIn)  # Assume full path missing if file doesn't exist?

        with open(fileIn, 'rb') as handle:
            dataIn = pickle.load(handle)

            if self.verbose['main']:
                print(f'Read data from {fileIn} with pickle.')

#         fInd = list(dataIn.keys())[-1]  # Set final key from data - probably not a robust method however! YES - see better method below
                                                  # This is currently used for indexing

        # Parse keys for max N - with checks in case list(dataIn.keys())[-1] fails (e.g. after data handling)
        fInd = np.array([k for k,v in dataIn.items() if isinstance(k,int)]).max()

        # Stack to data with new keys...?
        # Can use tuples, but might break other functions (which assume ints)?
        # YEP - currently breaks self.anayseFits()
#         self.data.update({((ind,k) if isinstance(k,int) else k):v for k,v in dataIn.items()})

        # Stack by ints only, but preserve other data per file
        self.data.update({((ind,k) if not isinstance(k,int) else k+fOffset):v for k,v in dataIn.items()})

        fOffset = fOffset + fInd + 1 # Update total N

        # If fits are batched, fix steps to batch size
        # This allows for missing cases in input file, and preserves batch step size
        # TODO: allow for batch size per file.
        if batch is not None:
            if np.mod(fOffset,batch):
                fOffset = fOffset + (batch - np.mod(fOffset,batch))

    self.fitInd = fOffset
