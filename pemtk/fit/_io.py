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
import pandas as pd

from epsproc import multiDimXrToPD, writeXarray, getFiles

# dataPath = Path(r'/home/jovyan/work/pemtk_fitting_runs_April2022/ioTests')
# dataPath = Path(r'/home/jovyan/ioTestsLocal/data_dumps')  # FF Docker

# import os

def setTimeStampedFileName(self,outStem = None, n = None, ext = 'pickle', timeString = None, timeFormat = "%d%m%y_%H-%M-%S"):
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

    timeFormat : Datatime format string, optional, default = "%d%m%y_%H-%M-%S"

    TODO: additional formatting options, data[key][item] naming option?

    """

    if outStem is None:
        outStem = 'PEMtk_data_dump'

    if timeString is None:
        timeString = dt.now()

    # Use n for multiple fit checkpointing - may want a more sophisticated routine for this however.
    if n is not None:
        fName = f'{outStem}_n{n}_{timeString.strftime(timeFormat)}'
    else:
        fName = f'{outStem}_{timeString.strftime(timeFormat)}'

    if ext is not None:
        fName = fName + f'.{ext}'

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
        fName = self.setTimeStampedFileName(outStem = outStem, n = n, ext = ext)


    fOut = Path(dataPath,fName)

    if fType == 'pickle':
        # pickleData(self,fOut)  # NOTE - here only include self for testing!
        sFlag = self._pickleData(fOut)

    elif fType == 'pdHDF':
        # _writePDData(self,fOut,**kwargs)  # NOTE - here only include self for testing!
        sFlag = self._writePDData(fOut,**kwargs)

    # elif fType == 'nc':
    # 30/06/22: updated handling on ePSproc.IO.writeXarray, so assume this as default case.
    else:
        # _writeXRData(self,fOut,**kwargs)  # NOTE - here only include self for testing!
        sFlag = self._writeXRData(fOut, engine=fType,**kwargs)

    if self.verbose['main']:
        if sFlag:
            print(f'Dumped data to {fOut} with {fType}.')
        else:
            print(f'Failed to write data to {fOut} with {fType}.')

    # Log file and return
    self.files['dataPaths'].append(dataPath)
    if sFlag:
        self.files['filesOut'].append(fOut)
    else:
        self.files['filesOutFailed'].append(fOut)

    return fOut


def processedToHDF5(self, dataKey = 'fits', dataTypes = ['dfLong','AFxr'], fType = 'pdHDF',
              outStem=None, multiFile = False, timeStamp = True, **kwargs):  #fOut, dataKey = 'fits', dataType='dfLong'):   # dataPath = None, fName = None, outStem = None):  # Assume fOut handled by wrapper
    """
    Save processed fit data to HDF5.

    Write self.data['fits']['dfLong'] and self.data['fits']['AFxr'] to file.

    Wraps self.writeFitData for processed data types.

    TODO: generalise to arb set of dataTypes and add checks.

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
        if not multiFile:
            outStemItem = outStem+dataType

        # Use outStem with or without timestamp (pass as fName)
        if timeStamp:
            self.writeFitData(dataKey = dataKey, dataType = dataType, outStem = outStemItem, fType = fType, **kwargs)

        else:
            self.writeFitData(dataKey = dataKey, dataType = dataType, fName = f'{outStemItem}.{fType}', fType = fType, **kwargs)


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

    return 1


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

    return 1



def _writeXRData(self, fOut, dataKey = 'fits', dataType='AFxr', **kwargs):
    """
    Dump item self.data[dataKey][dataType] to file using ep.writeXarray, built on Xarray.to_netcdf()

    This works well for basic Xarray structures.

    For complex data to netCDF need to either set engine='h5netcdf', forceComplex=True, or split to Re + Im groups (currently default for ep.writeXarray).
    For HDF5 complex data is supported.

    Attributes may also need to be sanitised for netCDF writer.

    See :py:func:`epsproc.IO.writeXarray`__ and https://docs.xarray.dev/en/latest/user-guide/io.html

    For docs on ePSproc Xarray IO (June 2022): https://epsproc.readthedocs.io/en/dev/dataStructures/ePSproc_dataStructures_IO_demo_280622.html

    Currently supported methods/engines:

    - netCDF via 'h5netcdf', 'scipy' or 'netcdf4' engines (using Xarray.to_netcdf())
    - 'hdf5' via h5py (via dictionary conversion routines, see :py:func:`epsproc.ioBackends.hdf5IO.writeXarrayToHDF5`__)

    """

    item = self.data[dataKey][dataType].copy()
    # item = item.drop('Euler')   # FOR TESTING ONLY!

    # writeXarray(item, fileName = fOut.name, filePath = fOut.parent, **kwargs)

    try:
        writeXarray(item, fileName = fOut.name, filePath = fOut.parent, **kwargs)
        return 1

    except Exception as e:
        print(f"\n*** Failed to write self.data[{dataKey}][{dataType}] to file {fOut} with ep.writeXarray().")
        print(f"ep.writeXarray() caught exception: {e}")
        print(f'(Partial write may have completed.)')
        return 0



#************ LOAD FILES

def getFilesList(fileIn = None, fileBase = None, fType = 'pickle', verbose = False):
    """Thin wrapper for epsproc.IO.getFiles - get file list from dir (no subdirs) by type"""

    return getFiles(fileIn = fileIn, fileBase = fileBase, fType = fType, verbose = verbose)


def loadFitData(self, fList = None, dataPath = None, batch = None, **kwargs):
    """
    Load data dumps from a file or set of files (and stack).

    Currently Pickle files only.

    See writeFitData for other options/file types - to be added here too.

    NOTE: currently only supports single dir for dir scan.
    For file lists can pass names only, in which case dataPath will be added, or pass full paths in list.

    """


    if dataPath is None:
        # dataPath = os.getcwd()  # OR
        dataPath = Path().absolute()

    if fList is None:
        # Run dir scan
        fList = getFilesList(fileBase = dataPath, fType = 'pickle')

    # NOTE: seems like bad logic here, but allows for full file paths in fList since ep.getFiles() further checks this!
    # Probably not very efficient however.
    else:
        if not isinstance(fList,list):
            fList = [fList]

        # Check passed list of files for validity
        fList = getFilesList(fileIn = fList, fileBase = dataPath, fType = 'pickle', verbose = self.verbose['sub'])

    if not fList:
        print('* Warning: empty file list for reader, please check paths.')

    # Check files exist - either on full path or cwd.
    # Already included in main loop below!
    # TODO: add error message here?
    # UPDATE: now also done above via getFilesList() functionality.


    fOffset = 0
    self.files['batches'] = {}
    for ind,fileIn in enumerate(fList):
    # for ind,fileIn in enumerate(fListChecked):
#         data.verbose['sub'] = 1

        # Now duplicated above, via getFilesList(), should remove here?
        # Except above returns strings, not Paths.
        if not Path(fileIn).exists():
            fileIn = Path(dataPath,fileIn)  # Assume full path missing if file doesn't exist?
        else:
            fileIn = Path(fileIn)  # Force Path object

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

        if batch is not None:
            # Stack by ints only, but preserve other data per file
            self.data.update({((ind,k) if not isinstance(k,int) else k+fOffset):v for k,v in dataIn.items()})
        else:
            # No batch case - just use inputs directly
            # self.data.update({k:v for k,v in dataIn.items()})
            self.data.update(dataIn)

        fStart = fOffset
        fOffset = fOffset + fInd + 1 # Update total N

        # If fits are batched, fix steps to batch size
        # This allows for missing cases in input file, and preserves batch step size
        # TODO: allow for batch size per file.
        if batch is not None:
            if np.mod(fOffset,batch):
                fOffset = fOffset + (batch - np.mod(fOffset,batch))

        # Log file input
        self.files['filesIn'].append(fileIn)
        self.files['dataPaths'].append(fileIn.parent)
        self.files['batches'][ind] = {'file': fileIn, 'fits':fInd, 'fitInds':[fStart, fOffset]}

    self.fitInd = fOffset
