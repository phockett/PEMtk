# PEMtk functions for data filtering
# Some similarities to existing ePSproc stuff, but specific to Pandas DataFrame for fit analysis (?)
# 28/09/21

# See ep.util.misc.checkDims - just added basic PD support.
# def checkDims(self, refDims)

from epsproc.util.misc import subselectDims

def thresFits(self, thres = None, dataType = None, key = 'fits', dataDict = 'dfLong'):
    """
    Very basic threshold mask for Pandas DataFrame. Note mask = True for values < thres.

    For more sophisticated multi-param filtering, see TMO-DEV filterData.
    Should be applicable here for Pandas data too...?
    """

    if dataType is not None:
        dimCheck = subselectDims(self.data[key][dataDict], refDims = dataType)

        mask = self.data[key][dataDict][dimCheck] < thres

    else:
        mask = self.data[key][dataDict] < thres

    self.data[key]['mask'] = mask


#************* FUNCTIONS FROM TMO-DEV
# NOT CURRENTLY IMPLEMENTED... but might want to use these for more flexibility.
# ALTHOUGH ALSO NOT GREAT/need some work!
# https://github.com/phockett/tmo-dev

def filterData(self, filterOptions = {}, keys = None, dim = 'energies', dTypes = ['raw','metrics']):
    """
    Very basic filter/mask generation function.

    filterOptions : dict, containing {dim:values} to filter on.
        Singular values are matched.
        Pairs of values as used as ranges.
        For multidim parameter sets, specify which source column to use as 3rd parameter.

    keys : list, optional, default = None
        Datasets to process, defaults to self.runs['proc']

    dim : str, optional, default = 'energies'
        Data to use as template. Not usually required, unless multidim return and/or default data is missing.

    dTypes : list, optional, default = ['raw','metrics']
        Data dicts to use for filtering.
        TODO: move this elsewhere!

    TODO:

    - More flexibility.
    - Filter functions, e.g. saturated electron detector shots? ('xc' > 0).sum(axis = 1) <=1000 in this case I think.

    07/12/20: added support for "metrics" data.

    """

    # Default to all datasets
    if keys is None:
        keys = self.runs['proc']

    for key in keys:

        # Check key/template dim exists
        # TODO: make this better.
        if dim not in self.data[key]['dims']:
            # Missing dim, try defaults.
            if 'gmd_energy' in self.data[key]['dims']:
                dim = 'gmd_energy'
            elif 'energies' in self.data[key]['dims']:
                dim = 'energies'

            # SACLA default fields for raw & process data (MIGHT CHANGE!)
            elif 'tagevent' in self.data[key]['dims']:
                dim = 'tagevent'

            elif 'tagID' in self.data[key]['dims']:
                dim = 'tagID'


            else:
                print(f"Missing dims for filtering: {dim}, and defaults not present.")

                # if self.verbose['main']:
                #     print(f"Reset )

        # Set full mask = true, use passed dim ('energies' as default) to define datasize event dim
        mask = np.ones_like(self.data[key]['raw'][dim]).astype(bool)
        if len(mask.shape)>1:
            mask = mask[:,0]

        for item in filterOptions.keys():

            # For 'raw' data types
            # if item in self.data[key]['raw'].keys():
            #     testData = np.array(self.data[key]['raw'][item])
            #
            # # For metrics (derived data)
            # elif item in self.metrics[key].keys():
            #     testData = self.metrics[key][item]

            # Version with dict testing
            # for dType in dTypes:
            #     if item in self.data[key][dType].keys():
            #         dataDict = dType
            #         testData = np.array(self.data[key][dataDict][item])  # np.array wrapper for 'raw' case
            #     else:
            #         dataDict = None
            #         testData = None

            testData = self.getDataDict(item, key, returnType = 'data')  # Method version of the above

            # if dataDict is not None:
            #     testData = np.array(self.data[key][dataDict][item])  # np.array wrapper for 'raw' case
            #
            # else:
            #     testData = None

                # 27/11/20 = quick hack for multidim col selection in v4 preprocessed data
                # Pass as dict with 'col', 'value' parameters
                # ACTUALY, multidim case below was OK, just had a bug!
                # if type(filterOptions[item]) is dict:
                #     col = filterOptions[item]['col']
                #     val = filterOptions[item]['val']
                #
                #     testData = testData[:,col]
                #

            if testData is not None:
                # Match single items
                if type(filterOptions[item]) is not list:
                    filterOptions[item] = [filterOptions[item]]  # UGLY - wrap to list.

                if len(filterOptions[item])==1:
                    mask *= (testData == filterOptions[item])

                if len(filterOptions[item])==2:
                    mask *= (testData >= filterOptions[item][0]) & (testData <= filterOptions[item][1])

                # Case for multidim testData
                if len(filterOptions[item])==3:
                    mask *= (testData[:,filterOptions[item][2]] >= filterOptions[item][0]) & (testData[:,filterOptions[item][2]] <= filterOptions[item][1])

            else:
                if self.verbose['main'] and (item is not 'desc'):   # Ignore 'desc' entries.
                    print(f"Can't filter on data type {item}, not found in dataset {key}")

            # TODO: add other filter types here.

        self.data[key]['mask'] = mask  # For single filter this is OK, for multiples see vmi version.


def getDataDict(self, dim, key = None, dTypes = None, returnType = 'dType', dropna = False):
    """
    Return specific dataset from various dictionaries by dimension name.

    dim : string
        Dimension (data) to find/check.

    key : string, int, optional, default = None
        Run key into main data structure.
        If None, use the first run in self.runs['proc'].

    dTypes : str, list, optional, default = self.dTypes
        Data dicts to check, defaults to global settings.

    returnType : str, optional, default = 'dType'
        - 'dType' return data type to use as index.
        - 'data' return data array.
        - 'lims' return min & max values only.
        - 'unique' return list of unique values.

    dropna : bool, optional, default = True
        Drop Nans in data?
        These cause issues for np.hist with 'auto' binning.
        But... in current code with basic masking, this breaks mask if array sizes are inconsistent.
        Better to filter out with ranges?

    08/12/20: first attempt, to replace repeated code in various base functions, and allow for multiple types (e.g. 'raw', 'metrics' etc.)

    TODO: may also want to add datatype to array conversion routine, since this will otherwise default to float64 and can be memory hungy.
    May also want to add chunking here too.

    TO FIX: dTypes checking buggy, for multiple matched dTypes only returns last matching item.

    """

    # Default to first dataset
    if key is None:
        key = self.runs['proc'][0]

    if dTypes is None:
        dTypes = self.dTypes

    dataDict = None # Set default
    for dType in dTypes:
        if (dType in self.data[key].keys()):  # && (dim in self.data[key][dType].keys()):  # Short circuit here in case dType doesn't exist in dict.
            if (dim in self.data[key][dType].keys()):
                dataDict = dType
            # else:
            #     dataDict = None

    if returnType is 'dType':
        return dataDict

    elif returnType is 'data':
        if dataDict is not None:
            if dataDict is 'raw':
                # return np.array(self.data[key][dataDict][dim])  # Explicit conversion to np.array - may just wrap everything this way?
                dataOut = np.array(self.data[key][dataDict][dim])
            else:
                # return self.data[key][dataDict][dim]
                dataOut = self.data[key][dataDict][dim]

            if dropna:
                return dataOut[~np.isnan(dataOut)]
            else:
                return dataOut

        else:
            return dataDict  # Currently set to None if dim not found, may change later.

    # def checkDims(testArray):
    #     """Check if array is 2D"""
