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

from epsproc import multiDimXrToPD, writeXarray

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


def writeFitData(self, dataPath = None, fName = None, outStem = None, n=None, fType = 'pickle', ext = None, **kwargs):
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

    elif fType == 'pdHDF':
        # _writePDData(self,fOut,**kwargs)  # NOTE - here only include self for testing!
        self._writePDData(fOut,**kwargs)

    elif fType == 'nc':
        # _writeXRData(self,fOut,**kwargs)  # NOTE - here only include self for testing!
        self._writeXRData(fOut,**kwargs)

#     if self.verbose['main']:
#         print(f'Dumped data to {fOut} with pickle.')

    return fOut


def aggToHDF5(self, dataKey = 'fits', dataTypes = ['dfLong','AFxr'], fType = 'pdHDF',
              outStem=None, multiFile = False, **kwargs):  #fOut, dataKey = 'fits', dataType='dfLong'):   # dataPath = None, fName = None, outStem = None):  # Assume fOut handled by wrapper
    """
    Save aggregate fit data to HDF5.

    Write self.data['fits']['dfLong'] and self.data['fits']['AFxr'] to file.

    Wraps self.writeFitData for processed data types.

    """

    if outStem is None:
        outStem = ''
    else:
        outStem = outStem + '_'

    if dataKey not in self.data.keys():
        print(f"*** Skipping data save, self.data[{dataKey}] missing.")

        return 0

    if 'dfLong' not in self.data[dataKey].keys():
        self.analyseFits()

    # Save dfLong
    # self.data['fits']['dfLong'].to_
    # self.writeFitData(self, dataPath = dataPath, outStem = 'pd_dump_test', fType = 'pdHDF')

    # Write to file using writeFitData wrapper.
    # Include outStem=dataType to avoid accidental datastamp file overwrite, or could support passed outStem?
    # UPDATE - now implemented above, set multiFile=True to use
    for dataType in dataTypes:
        # self.writeFitData(self, dataKey = dataKey, dataType = dataType, fType = fType, **kwargs)
        outStemItem = outStem
        if multiFile:
            outStemItem = outStem+dataType

        self.writeFitData(dataKey = dataKey, dataType = dataType, outStem = outStemItem, fType = fType, **kwargs)


#************* FILE WRITERS

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

def _writePDData(self, fOut, dataKey = 'fits', dataType='dfLong'):
    """
    Dump item self.data[dataKey][dataType] to file using Pandas.to_hdf().

    This works well for Pandas Dataframes, including complex data.
    Also works for Xarray (via Pandas conversion routine), but may lose attribs.

    Default case will append to file if it exists, so multiple calls will add items (nodes/groups) to hdf store.
    (Note read in via Pandas requires per key reading in this case.)

    See https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_hdf.html

    And https://pandas.pydata.org/docs/user_guide/io.html#hdf5-pytables

    """

    item = self.data[dataKey][dataType]

    if not isinstance(item, pd.DataFrame):
        try:
            item, _ = multiDimXrToPD(item, colDims = 'Fit', squeeze = False)
        except:
            print(f"*** Couldn't convert self.data[{dataKey}][{dataType}] to Pandas, skipping file writing")
            return 0

    item.to_hdf(fOut, dataType)   # Write to file with named group dataType

    if self.verbose['main']:
        print(f'Dumped self.data[{dataKey}][{dataType}] to {fOut} with Pandas .to_hdf() routine.')



def _writeXRData(self, fOut, dataKey = 'fits', dataType='AFxr', **kwargs):
    """
    Dump item self.data[dataKey][dataType] to file using ep.writeXarray, built on Xarray.to_netcdf()

    This works well for basic Xarray structures.

    For complex data need to either set engine='h5netcdf', forceComplex=True (currently default for ep.writeXarray), or split to Re + Im groups.

    Attributes may also need to be sanitised for netCDF writer.

    See :py:func:`epsproc.IO.writeXarray` and https://docs.xarray.dev/en/latest/user-guide/io.html


    """

    item = self.data[dataKey][dataType].copy()
    # item = item.drop('Euler')   # FOR TESTING ONLY!

    # writeXarray(item, fileName = fOut.name, filePath = fOut.parent, **kwargs)

    try:
        writeXarray(item, fileName = fOut.name, filePath = fOut.parent, **kwargs)

    except Exception as e:
        print(f"*** Failed to write self.data[{dataKey}][{dataType}] to file with ep.writeXarray().")
        print(e)



#************ LOAD FILES


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
