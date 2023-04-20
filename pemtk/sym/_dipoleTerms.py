"""
Compute dipole terms and allowed product functions for libmsym and symmetrized harmonics routines.

17/04/23 v1

"""
from copy import deepcopy
import xarray as xr
import numpy as np

from epsproc import multiDimXrToPD, matEleSelector
from ._util import toePSproc

# Generate basis

# V1 STAND-ALONE VERSION
# Note this uses symHarm class
#
# from pemtk.sym.symHarm import symHarm
#
# def dipoleTerms(PG='Cs', reCartMap = {-1:'y',0:'z',1:'x'}, dimMap = {'C':'Cont', 'mu':'it'}):
#     """Compute dipole terms from symHarm object."""
#
#     # Create object
#     sym = PG
#     dipSym = symHarm(PG = sym, llist=1)  # In future can use llist = 1 to compute only l=1 term!
#
#
#     # Set (x,y,z) labels from real harmonics
#     dfCart = dipSym.coeffs['DF']['real'].rename(index = reCartMap).droplevel(['l','h','mu'])
#     cartMapping = dfCart.index.to_list()
#
#     # Set m terms from comp harmonics
#     # Select object and process - no remappinng
#     # dipolePD, _ = multiDimXrToPD(dipSym.coeffs['XR']['b (comp)'], colDims = 'C', thres=1e-4) #, fillna=True)
#
#     # With remapping - optional?
#     dataType = 'matE'
#     dipSym.toePSproc(dimMap = dimMap, dataType=dataType)
#     dipolePD, _ = multiDimXrToPD(dipSym.coeffs[dataType]['b (comp)'], colDims = 'Cont', thres=1e-4) #, fillna=True)
#
#     dictOut = {}
#     # for (col, vals) in dipSym.coeffs['DF']['comp'].copy().unstack(level='C').iteritems():
#     for (col, vals) in dipolePD.iteritems():
#         dictOut[col]={'m':vals.dropna().index.get_level_values('m').to_list(),
#                       'pol':[p[1] for p in cartMapping if p[0]==col]}
#
#     return dipSym, dictOut


# TODO: generalise this to any pair of symmetries. Then chain for allowed terms for orbital etc.
# MAY want to do this as (1) multiple table look-ups, (2) per directProductTable() for specific syms only (mulitple times) - might be cleaner?
# def allowedProducts(dipSym, dictOut):
#     """Compute allowed continuum X dipole terms."""
#
#     # Check direct products
#     # dfDP,_ = diretProductTable(PG = dipSym.PG)
#     dfDP = dipSym.directProducts
#
#     prodTerms = dfDP[list(dictOut.keys())]  # Subselect valid cols
#
#     # Check which terms contain totally symmetric representation
#     # prodTerms.applymap(lambda x: True if (dipSym.ctx.character_table.symmetry_species[0].name in x) else False)
#     prodValid = prodTerms.applymap(lambda x: True if (dipSym.irreps[0] in x) else False)
#
#     # List of tuples of valid (C,dip) pairs
#     dipContList = prodValid[prodValid > 0 ].stack().index.tolist()
#
#     return dipContList, prodTerms, prodValid

# v2 Similar to above, but for symHarm class - this assumes self and subselects on l.
def dipoleTermsSymHarm(self, reCartMap = {-1:'y',0:'z',1:'x'}, dimMap = {'C':'Cont', 'mu':'it'}):
    """Compute dipole terms from symHarm object."""

    # Create object
    # sym = PG
    # dipSym = symHarm(PG = sym, llist=1)  # In future can use llist = 1 to compute only l=1 term!

    # dipSym = deepcopy(self)  # For class version set copy to allow for subselection etc. below.
                               # Can't deepcopy class, so just set relevant parts below instead.

    # Set (x,y,z) labels from real harmonics, l=1
    dfCart = self.coeffs['DF']['real'].xs(1, level='l').rename(index = reCartMap).droplevel(['h','mu'])
    cartMapping = dfCart.index.to_list()

    # Set m terms from comp harmonics
    # Select object and process - no remappinng
    # dipolePD, _ = multiDimXrToPD(dipSym.coeffs['XR']['b (comp)'], colDims = 'C', thres=1e-4) #, fillna=True)

    # With remapping - optional?
    dataType = 'matE'
    # dipSym.toePSproc(dimMap = dimMap, dataType=dataType)
    coeffs= toePSproc(self.coeffs, dimMap = dimMap, dataType = dataType, verbose = self.verbose)
    dipoleTerms = matEleSelector(coeffs['b (comp)'].unstack(), inds = {'l':1}, thres = None, drop=True)
    dipoleTerms.attrs = coeffs.attrs  # Propagate attribs
    dipoleTerms.attrs['dataType'] = 'dipole'

    dipolePD, _ = multiDimXrToPD(dipoleTerms, colDims = 'Cont', thres=1e-4) #, fillna=True)

    dictOut = {}
    # for (col, vals) in dipSym.coeffs['DF']['comp'].copy().unstack(level='C').iteritems():
    for (col, vals) in dipolePD.iteritems():
        dictOut[col]={'m':vals.dropna().index.get_level_values('m').to_list(),
                      'pol':[p[1] for p in cartMapping if p[0]==col]}

    # return dipSym, dictOut
    # return {'dictOut':dictOut, 'dipolePD':dipolePD, 'dipoleTerms':dipoleTerms, 'coeffs':coeffs, }
    self.dipole = {'dipoleSyms':dictOut, 'dipolePD':dipolePD, 'dipoleXR':dipoleTerms, 'coeffsXR':coeffs, }

    if self.verbose:
        print("Found dipole symmetries: ")
        print(dictOut)


# As above, but use self.
def allowedProductsTable(self):
    """Compute allowed continuum X dipole terms."""

    # Check direct products
    # dfDP,_ = diretProductTable(PG = dipSym.PG)
    dfDP = self.directProductTable

    prodTerms = dfDP[list(self.dipole['dipoleSyms'].keys())]  # Subselect valid cols

    # Check which terms contain totally symmetric representation
    # prodTerms.applymap(lambda x: True if (dipSym.ctx.character_table.symmetry_species[0].name in x) else False)
    prodValid = prodTerms.applymap(lambda x: True if (self.irreps[0] in x) else False)

    # Drop other terms?
    prodValid = prodValid.where(prodValid).dropna(how='all').fillna('')

    # List of tuples of valid (C,dip) pairs
    dipContList = prodValid[prodValid == True ].stack().index.tolist()

    # return dipContList, prodTerms, prodValid
    self.dipole.update({'dipContList':dipContList, 'prodTerms':prodTerms, 'prodValid':prodValid})


# Assign allowed subset and m values via symmetry
def assignSymMuTerms(self, keyDim = 'Cont',    # targSym = None,
                    dataType = 'matE', dimMap = None,
                    dataTypeOut = 'symAllowed'):
    """
    Assign matrix elements from dipole-allowed symmetries, and associated mu values.

    20/04/23 v2 - use self.continuum instead of basic direct product for assignments
    18/04/23 v1

    """

    # Convert to ePSproc data type if not already present
    if not dataType in self.coeffs.keys():
        if dimMap is None:
            self.toePSproc(dataType=dataType)
        else:
            self.toePSproc(dataType=dataType, dimMap = dimMap)  # With custom dimMap

    # Set to totally symmetric target if not specified
    # if targSym is None:
    #     targSym = self.irreps[0]

    # Set m from keyDim and allowed terms...
    testAllowed = {}
    dipoleTerms = self.dipole['dipoleSyms'].copy()
    mMapping = {}

    symTerms = self.coeffs[dataType].coords[keyDim].values

# v1
#     # Compute as direct products
#     # Basically already in self.dipole['prodValid']... so could use that instead...
#     for n, s1 in enumerate(symTerms):
#     # for n, s1 in enumerate(testDP):
#         for s,v in dipoleTerms.items():
#             testAllowed[(s1,s)] = {'directProd':self.directProduct(terms = [s1, s])}
#             testAllowed[(s1,s)]['allowed'] = True if self.irreps[0] in testAllowed[(s1,s)]['directProd'] else False

#             if testAllowed[(s1,s)]['allowed']:
#                 mMapping[s1] = v

#********** v2 - use self.continuum
    testAllowed = self.continuum['dict']

    for k,v in testAllowed.items():
        # print(v)
        if v['allowed']:
            mMapping[v['targSym']] = v

#**********

    # Set mu axis correctly...
    muList = []
    testDS = self.coeffs[dataType].copy()
    for k,v in mMapping.items():
        muXR = testDS.where(testDS.coords[keyDim] == k, drop=True)

        if len(v['m']) > 1:
            muXR0 = muXR.copy()
            muXR0.coords['mu'] = muXR.mu.fillna(v['m'][0])
            muXR1 = muXR.copy()
            muXR1.coords['mu'] = muXR.mu.fillna(v['m'][1])

            muList.extend([muXR0.copy(),muXR1.copy()])

        else:
            muXR.coords['mu'] = muXR.mu.fillna(v['m'])
            muList.extend([muXR.copy()])

    testMuMerge = xr.merge(muList)  # THIS MIGHT BE OK, not sure yet...
    mergePD, _ = multiDimXrToPD(testMuMerge['b (comp)'], colDims = keyDim, thres=1e-4)

    self.coeffs[dataTypeOut] = {'XR':testMuMerge, 'PD':mergePD}

    if self.verbose:
        print(f"Assigned dipole-allowed terms for dim = '{keyDim}' to self.coeffs['{dataTypeOut}']")


def assignMissingSym(self, dim, values, dataType = 'matE', multiInd = 'Sym'):
    """
    Replace values in dim, for levels of a (Symmetry style) MultiIndex.

    Use this to change symmetry labels, e.g. to replace missing symmetries assigned as 'U'.

    Parameters
    ----------
    dim : str
        Dimension to replace, from matE dataType.
        Can be ['Cont','Targ','Total']

    values : str, list or array
        If string, all values are replace by this.
        If list or array must match size of existing coords.

    dataType : str, optional, default = 'matE'
        DataType to use, from self.coeffs[dataType]
        Must correspond to Xarray variable with multiInd present.

    multiInd : str, optional, default = 'Sym'
        MultiIndex containing levels to replace.
    """

    newDS = self.coeffs[dataType].copy()

    if isinstance(values, str):
        newVals = newDS.coords[dim].str.replace('U',values).values   # OK
    else:
        newVals = values


    # Quick XARRAY MultiIndex COORD REFRESHER
    # see https://github.com/pydata/xarray/discussions/6936?sort=top
    # Think I've done this elsewhere - need to unstack and restack for multiindex coord replacement.

    # testDS = testDS.reset_index(dim, drop=True).assign_coords(Targ=('Sym',newTarg.values))   # .set_index(Sym='Targ', append=True)   # Not required?
    newDS = newDS.reset_index(dim, drop=True).assign_coords({dim:(multiInd,newVals)}).set_index(Sym=dim, append=True)  # Output looks OK without set_index too, but seems to cause issues later with dim dropping!

    self.coeffs[dataType] = newDS

    if self.verbose:
        print(f"*** Updated self.coeffs['{dataType}'] with new coords.")


def assignMissingSymProd(self, dim = 'Total', dataType = 'matE', multiIndName = 'Sym'):
    """
    Assign missing symmetry label as direct product of other existing labels.

    For a given dimension, assign values from other dims of the same multiindex set.

    Default case uses dim = 'Total', dataType = 'matE', multiIndName = 'Sym'.
    Sym contains dims 'Sym': ['Cont', 'Targ', 'Total']
    (See ep.dataTypesList())

    """

    # Check dims?
    # ep.util.misc.checkDims
    # ep.dataTypesList()

    # Assume dims...?
    # dims = ep.util.listFuncs.getRefDims(dipSym.coeffs[dataType])
    dims = set(self.coeffs[dataType].indexes[multiIndName].names) - {dim}  # Get dims from set

    # Case for symAllowed only
    # cont = dipSym.coeffs['symAllowed']['XR'].indexes['Sym'].get_level_values('Cont')
    # targ = dipSym.coeffs['symAllowed']['XR'].indexes['Sym'].get_level_values('Targ')

    # Assign all cases
    # This is OK, but may want to iterate over each set instead to allow >2 later...?
    symsInput = {}
    symsList = []
    for d in dims:
        symsInput[d] = self.coeffs[dataType].indexes[multiIndName].get_level_values(d)

        symsList.append(list(symsInput[d]))

    # Use np.array to handle arb sets and transpose
    prodList = np.array(symsList).T

#     symsInput = {}
#     symsList = []

#     d1 = self.coeffs[dataType].indexes[multiIndName].get_level_values(dims[0])
#     for n,item in enumerate(d1):
#         for d in dims:
#             symsInput[d] = self.coeffs[dataType].indexes[multiIndName].get_level_values(d)

    # targ = dipSym.coeffs[dataType].indexes['Sym'].get_level_values('Targ')


    # for item in dipSym.coeffs['symAllowed']['XR'].indexes['Sym']:
    #     # print(item[1])
    #     newList.append([*item]).str.replace('U','test')
    #     # newList.append(


    # return symsInput, symsList

    directProdDict = {}
    directProdList = []

    for n,item in enumerate(prodList):
        directProd = self.directProduct(terms = item)
        # newDict[n] = {'cont':cont[n], 'targ':targ[n], 'prod':directProd}
        directProdDict[n] = {'Terms':item, 'Product':directProd}
        directProdList.append(directProd[0])

        if len(directProd) > 1:
            print(f"Found multiple products for {' x '.join(item)} = {directProd}; assigning '{dim}' as {directProd[0]}.")
        else:
            print(f"Assigned '{dim}' from {' x '.join(item)} = {directProd}")


    self.assignMissingSym(dim,directProdList)  # Assign dim with list
