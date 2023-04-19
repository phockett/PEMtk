"""
Direct product functions for libmsym and symmetrized harmonics routines.

17/04/23 v1

Developed/adapted from libmsym demo code, see https://github.com/mcodev31/libmsym/issues/29
Special thanks to mcodev31 for the help.

For quick reference tables, see also http://www.gernot-katzers-spice-pages.com/character_tables/index.html

Formulas should follow those on p29 in Atkins, B.P.W., Child, M.S. and Phillips, C.S.G. (2006) Tables for Group Theory. Oxford University Press.

"""

import libmsym as msym, numpy as np, pandas as pd

#********* LOW LEVEL FUNCTIONS
def prod(character_table,s1,s2):
    """
    Array product of two irrep characters.

    Parameters
    ----------
    character_table : libmsym character table.

    s1 : int
        Index for irrep/symmetry species 1.

    s2 : int
        Index for irrep/symmetry species 2.

    Returns
    -------
    np.array : array product.

    """

    # Looped style
    # r = np.zeros(len(character_table.symmetry_operations));
    # for i in range(len(character_table.symmetry_operations)):
    #     r[i] = character_table.table[s1][i]*character_table.table[s2][i]

    # Numpy array style - just need array product of representations here
    r = character_table.table[s1]*character_table.table[s2]
    return r;


def prodList(character_table,symList):
    """
    Array product of (any number of) irrep characters from a list.

    Parameters
    ----------
    character_table : libmsym character table.

    symList : list
        List of species to use.
        E.g. symList = ['A1g','A1u','B1u', 'B1u'] defines resultant as direct product of these terms.

    Returns
    -------
    np.array : array product.


    Examples
    --------
    >>> PG = 'D2h'
    >>> ctx = msym.Context(elements=[msym.Element(name = "H", coordinates = [.0,.0,.0])], basis_functions=[], point_group=PG)
    >>> r = prodMulti(ctx.character_table, ['A1g','A1u','B1u', 'B1u'])
    >>> r
    array([ 1.,  1., -1., -1., -1., -1.,  1.,  1.])

    """

    # Symmetry lookup table.
    symLookup = {item.name:n for n,item in enumerate(character_table.symmetry_species)}

    # Numpy array style - just need array product of representations here
    # UPDATE: now assigned in loop.
    # r = np.ones(len(character_table.table[0]))
    # r = character_table.table[0]

    for n,s in enumerate(symList):
        # TODO: add try/except for error checking here.
        sInd = symLookup[s]

        if n == 0:
            r = character_table.table[sInd].copy()
        else:
            r *= character_table.table[sInd]

    return r;


# Original version (looped)
# def dot(character_table,s,p):
#     """Irrep representation dot products.""""
#
#     r = 0
#     l = 0
#     for i in range(len(character_table.symmetry_operations)):
#         l += character_table.class_count[i];
#         r += p[i]*character_table.table[s][i]*character_table.class_count[i]
#     return int(round(r))//l;

# Rewrite array style
def dotArray(character_table,s,p):
    """Irrep representation dot products."""
    h = np.array(character_table.class_count)
    r = p * character_table.table[s] * h

    return int(round(sum(r)))//sum(h)
    # return int(round(r))//h.sum();
    # return r, h

# Try the final analysis/projection/dot products independently.
def symDecomposition(character_table, p):
    """
    Decompose basis vector to irreps

    Parameters
    ----------
    character_table : libmsym character table.

    p : np.array
        Vecotor for analysis.
        Must contain len(character_table.symmetry_species) elements.

    """
    # Try array/dict style...
    dDict = {}
    # dArr = []  # List
    dArr = np.zeros(len(character_table.symmetry_species));  # Array

    for (k,s3) in enumerate(character_table.symmetry_species):
    #     # d = dot(ctx.character_table,k,p)

        dDict[s3.name] = dotArray(character_table,k,p)   # This will output full vector with labels
        # dArr.append(dotArray(ctx.character_table,k,p))   # Vector only
        # dArr[k] = dotArray(ctx.character_table,k,p)
        dArr[k] = dDict[s3.name].copy()

    return dDict, dArr


#******** TOP LEVEL FUNCTIONS

def diretProductFromList(PG="Cs", terms=[]):
    """
    Compute direct products from a list of irreps/species with [libmsym](https://github.com/mcodev31/libmsym).

    Developed from https://github.com/mcodev31/libmsym/issues/29.

    Parameters
    ----------
    PG : point group string.
        Supported point groups from libmsym: Ci, Cs, Cnv, Dn, Dnh, Dnd, Td, O, Oh, I and Ih

    terms : list
        List of species to use.
        E.g. symList = ['A1g','A1u','B1u', 'B1u'] defines resultant as direct product of these terms.

    Returns
    -------
    list : resultant as list of specices.

    dict : contains various other forms of the output (as labelled).


    Examples
    --------
    >>> dDict, _ = diretProductFromList(PG = 'D2h', terms = ['A1g','A1u','B1u', 'B1u'])
    >>> dDict
    {'A1g': 0,
     'B1g': 0,
     'A1u': 1,
     'B1u': 0,
     'B3g': 0,
     'B3u': 0,
     'B2g': 0,
     'B2u': 0}

    """

    # Create libmsym object
    ctx = msym.Context(elements=[msym.Element(name = "H", coordinates = [.0,.0,.0])], basis_functions=[], point_group=PG)

    # Vector product of specified terms
    p =  prodList(ctx.character_table, terms)

    # Decompose to irreps
    dDict, dArr = symDecomposition(ctx.character_table, p)
    resultsDict = {k:v for k,v in dDict.items() if v > 0}

    return list(resultsDict.keys()), {'resultsDict':resultsDict,
            'fullDict':dDict,'fullArr':dArr, 'vectorProduct':p, 'libmsymObj':ctx}


def directProductTable(PG="Cs"):
    """
    Generate direct product tables with [libmsym](https://github.com/mcodev31/libmsym).

    Developed from https://github.com/mcodev31/libmsym/issues/29


    Parameters
    ----------
    PG : point group string.
        Supported point groups from libmsym: Ci, Cs, Cnv, Dn, Dnh, Dnd, Td, O, Oh, I and Ih

    Returns
    -------
    pd.DataFrame : direct product table.

    dict : contains various other forms of the output (as labelled).

    """
    dpDict = {}
    tabOut = []
    symList = []

    # with msym.Context(elements=[msym.Element(name = "H", coordinates = [.0,.0,.0])], basis_functions=[], point_group="D2h") as ctx:
    # Remove context manager to pull out ctx
    ctx = msym.Context(elements=[msym.Element(name = "H", coordinates = [.0,.0,.0])], basis_functions=[], point_group=PG)
    symLabels = [item.name for item in ctx.character_table.symmetry_species]
    for (i,s1) in enumerate(ctx.character_table.symmetry_species):
        rowVec = []
        rowTuple = []
        rowSimple = []
        for (j,s2) in enumerate(ctx.character_table.symmetry_species):
            p = prod(ctx.character_table,i,j);
            # print(p)
            # print(s1.name+'x'+s2.name+'=',end='')
            # s = ''
            # for (k,s3) in enumerate(ctx.character_table.symmetry_species):
            #     # d = dot(ctx.character_table,k,p)
            #     d = dotArray(ctx.character_table,k,p)
            #     # print(d)
            #     # print(k)
            #     if(d != 0):
            #         if s != '':
            #             s += '+'
            #         s += str(d) + s3.name

            # Try array/dict style...
            dDict = {}
            # dArr = []  # List
            dArr = np.zeros(len(ctx.character_table.symmetry_species));  # Array
            for (k,s3) in enumerate(ctx.character_table.symmetry_species):
            #     # d = dot(ctx.character_table,k,p)

                dDict[s3.name] = dotArray(ctx.character_table,k,p)   # This will output full vector with labels
                # dArr.append(dotArray(ctx.character_table,k,p))   # Vector only
                dArr[k] = dotArray(ctx.character_table,k,p)

            # print(dDict)
            # print(dArr)
            # print(s)

            # dpDict.update({(s1.name, s2.name):s})   # Use multiindex style
            rowVec.append(dArr)  # All vectors
            rowTuple.append([(m, symLabels[n]) for n,m in enumerate(dArr) if m>0])  # Keep only non-zero values, format as tuples
            rowSimple.append([symLabels[n] for n,m in enumerate(dArr) if m>0])  # List only allowed terms


        # dpDict.update({s1.name:row})  # Row/col style
        dpDict.update({s1.name:rowTuple})

        # List style
        # tabOut.append(row)
        tabOut.append(rowSimple)
        symList.append(s1.name)

    # Create DataFrame from results
    # df = pd.DataFrame(tabOut, index = symList, columns = symList)  # Tabulate as full output
    df = pd.DataFrame(tabOut, index = symList, columns = symList)  # Tabulate as table

    return df, {'ctx':ctx,'rowVec':rowVec,'rowTuple':rowTuple,'rowSimple':rowSimple,'dpDict':dpDict,'tabOut':tabOut,'symList':symList}
