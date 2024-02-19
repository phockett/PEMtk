"""
Class for determining and handling symmetrized harmonics.

23/03/22    v2

- Switched to nested dicts for coeffs, now all in self.coeffs[dtype].

24/02/22    v1 in development.

- For early dev work see http://localhost:8888/lab/workspaces/symm/tree/python/chem/tools/symmetrized_harmonics_libmsym_tests_160122.ipynb
- For class dev see http://localhost:8888/lab/workspaces/symm/tree/python/chem/tools/symmetrized_harmonics_PEMtk-dev_240222.ipynb

TODO

- Data to dicts?  Currently have multiple self.coeffXX attribs.
- Interface/convert to matrix elements for ePSproc use.

"""

import libmsym as msym

import numpy as np
import pandas as pd
import argparse, random, sys

import pyshtools as pysh

try:
    import xarray as xr
    xrFlag = True

except:
    print("***Xarray not found, XR outputs not available.")
    xrFlag = False

from ._util import prefixAxisTitles, listPGs, toePSproc
from ._directProduct import directProductTable, diretProductFromList
# from ._dipoleTerms import dipoleTerms, allowedProducts


class symHarm():
    r"""
    Class to set and manipulate symmetrized harmonics.

    Main symmetry routine uses [libmsym](https://github.com/mcodev31/libmsym), code adapted from https://github.com/mcodev31/libmsym/issues/21


    Background
    ----------

    Symmetrized (or generalised) harmonics, which essentially provide correctly symmetrized functions for a given
    irreducible representation, $\Gamma$, be defined by linear combinations of spherical harmonics \cite{Altmann1963a,Altmann1965,Chandra1987}:

    \begin{equation}
    X_{hl}^{\Gamma\mu*}(\theta,\phi)=\sum_{\lambda}b_{hl\lambda}^{\Gamma\mu}Y_{l,\lambda}(\theta,\phi)\label{eq:symm-harmonics}
    \end{equation}

    where:

    - $\Gamma$ is an irreducible representation,
    - ($l$, $\lambda$) define the usual spherical harmonic indicies (rank, order)
    - $b_{hl\lambda}^{\Gamma\mu}$ are symmetrization coefficients,
    - index $\mu$ allows for indexing of degenerate components,
    - $h$ indexe cases where multiple components are required with all other quantum numbers identical.

    The exact form of these coefficients will depend on the point-group of the system, see, e.g. \cite{Chandra1987,Reid1994}.


    Method
    ------

    - Point-groups, character table generation and symmetrization (computing $b_{hl\lambda}^{\Gamma\mu}$ parameters) is handled by [libmsym](https://github.com/mcodev31/libmsym) (based on https://github.com/mcodev31/libmsym/issues/21).
        - Supported point groups: Ci, Cs, Cnv, Dn, Dnh, Dnd, Td, O, Oh, I and Ih
        - Routines have been tested with libmsym v0.2.2, March 2022 (last commit https://github.com/mcodev31/libmsym/commit/c99470376270db4ec4c925b952fa722e011377d6).

    - Spherical harmonic expansions and conversions (real, imaginary, normalization etc.) and basic ploting (2D maps) are handled by [pySHtools](https://shtools.oca.eu).
        - Routines have been tested with v4.5 and 4.9 (current March 2022).

    - TODO: Interface to PEMtk/ePSproc for other plotters and handling routines.
    - TODO: some hardcoded labelling to fix, also incorrect in places ("character" instead of "irrep")


    Examples
    --------

    >>> # Compute params for Td, lmax=4
    >>> xlm = symHarm('Td',4)

    >>> # Print charater table
    >>> xlm.printCharacterTable()

    >>> # Display table of params
    >>> xlm.displayXlm()

    >>> # Plot functions
    >>> xlm.plotXlm()


    """

    from ._dipoleTerms import dipoleTermsSymHarm as dipoleTerms
    from ._dipoleTerms import allowedProductsTable, assignSymMuTerms, assignMissingSym, assignMissingSymProd

    def __init__(self, PG = 'Cs', lmax = 4, llist = None, dims = ['C', 'h', 'mu', 'l', 'm']):
        """
        Init class object.

        Paramters
        ---------
        PG : string, default = 'Cs'
            Point-group.

        lmax : int, default = 4
            Maximum l to use.

        llist : list, optional, default = None
            If set, used instead of lmax to define specific l to use (all m).
            E.g. llist = [1,3,5].
            TODO: more options here, e.g. set m, pass existing Xarray or SHtools objects?
            Can also subselect from outputs of course.

        dims : list, default = ['C', 'h', 'mu', 'l', 'm']
            Dimension labels to use for outputs.

        """

        # Set params
        self.PG = PG
        self.lmax = lmax
        self.llist = llist

        self.dims = dims # Dims names from defaults, or as passed.
        self.dtypes = ['O',int,int,int,int]  # Set dtypes for dims

        self.verbose = 1
        self.PGlist = listPGs()  # Set via function to allow for access from other classes etc.

        # Set data structure
        self.coeffs = {}

        # Set single element as expansion centre
        self.elements = [msym.Element(name = "C", coordinates = [0.0, 0.0,0.0])]

        # Compute harmonic coeffs
        self.set_lm_basis()
        self.calcSymHarmonics()

        # Conversions (with defaults)
        self.setCoeffsSH()
        self.setCoeffsXR()

        # Additional outputs
        self.getIrreps()
        self.directProductTable, self.directProductsDict = directProductTable(PG = self.PG)
        # self.directProduct = diretProductFromList

        # Dipole terms
        # self.dipoleTerms = dipoleTerms
        self.dipoleTerms()
        self.allowedProductsTable()
        # self.dipoleTerms = dipoleTerms(self.PG)  # TODO: may want to set dimMap here too?
        # self.dipoleProducts = allowedProducts(*self.dipoleTerms)



    #*********************** CALCULATION FUNCTIONS

    def set_lm_basis(self):
        """Set basis functions for each (element, l, m)."""

        self.basis_functions = []

        if self.llist is not None:
            # Wrap int case.
            if not isinstance(self.llist, list):
                self.llist = [self.llist]

            llist = self.llist

        else:
            llist = range(0, self.lmax+1)

        for l in llist:
            for m in range(-l,l+1):
        #         basis_functions.extend([set_basis(e, l=l, m=m, n=n, name=f"{l}s{m}") for e in elements])
                self.basis_functions.extend([self.set_basis(e, l=l, m=m, n=self.lmax+1, name=f"{l},{m}") for e in self.elements])


    # Set basis for arb l,m,n...
    def set_basis(self,element, l, m, n, name):
        """Set basis function for a given (element, n, l, m,)."""
        basis_function = msym.RealSphericalHarmonic(element = element, l=l, m=m, n=n, name = name)
        element.basis_functions = [basis_function]
        return basis_function


    def calcSymHarmonics(self):
        """Compute symmetrized harmonic coeffs for point group up to lmax."""

        # Try tabulting for Pandas
        tabOut = []

        ctx = msym.Context(elements = self.elements, basis_functions = self.basis_functions, point_group = self.PG)
        for ss in ctx.subrepresentation_spaces:
        #     print(ctx.character_table.symmetry_species[ss.symmetry_species].name)
            for salcix, salc in enumerate(ss.salcs):
        #         print(f"salc {salcix}")
                for pfix, pf in enumerate(salc.partner_functions):
        #             print(f"\tpartern function {pfix}")
                    for ix, v in enumerate(pf):
                        if (abs(v) > 100*sys.float_info.epsilon):
        #                     print("\t\t{0}*{1} ".format(v,salc.basis_functions[ix].name));
                            l,m = salc.basis_functions[ix].name.split(',')

                            # Set line in output - lists
                            # tabOut.append([ctx.character_table.symmetry_species[ss.symmetry_species].name, salcix, pfix, int(l), int(m), v])
                            # Tuple version for NP structured array creation.
                            tabOut.append((ctx.character_table.symmetry_species[ss.symmetry_species].name, salcix, pfix, int(l), int(m), v))

        # Log outputs - may want to use dicts?
        self.ctx = ctx
        # self.coeffTable = tabOut
        # self.coeffs['libmsym'] = {'real':tabOut}

        # Convert to NP structured array (keeps heterogeneous dtypes)
        dimTypes = list(zip(self.dims,self.dtypes))  # Dims
        dimTypes.extend([('b',float)])   # Data
        self.coeffs['libmsym'] = {'real':np.asarray(tabOut, dtype=dimTypes)}

        # Set PD table form too
        self.setCoeffsPD()


    def directProduct(self, terms):
        """
        Compute direct products for a list of terms (symmetries/species/irreps).
        """

        prodList, *_ = diretProductFromList(PG = self.PG, terms = terms)

        return prodList


    def directProductDipole(self, terms):
        """
        Compute direct products for a list of terms (symmetries/species/irreps), and include dipole terms.
        """

        # Dipole terms
        dipoleDict = self.dipole['dipoleSyms']

        # Loop over dipole terms and find direct products
        dipDictOut = {}
        for s in dipoleDict:
            termsD = terms.copy()
            termsD.append(s)
            # print(termsD)
            rList, *_ = diretProductFromList(PG = self.PG, terms=termsD)
            # dipDict[s] = rDict.copy()
            dipDictOut[s] = dipoleDict[s].copy()
            dipDictOut[s]['dipSym'] = s
            dipDictOut[s]['result'] = rList.copy()
            # dipDict[s]['valid'] = [k for k in rDict.keys() if k == 'A1g'] # self.irreps[0]]
            dipDictOut[s]['allowed'] = True if self.irreps[0] in rList else False  # Allowed terms contain totally symmetric rep.
            dipDictOut[s]['terms'] = terms

        self.dipole['dipoleProducts'] = dipDictOut.copy()


    def directProductContinuum(self, terms, disp = True):
        """
        Compute direct products for a list of terms (symmetries/species/irreps), and include dipole terms.

        Similar to directProductDipole, but additionally loops over all irreps to determine valid photoionization combinations, i.e. allowed continuum symmetries for the given cases.
        """

        # Dipole terms
        dipoleDict = self.dipole['dipoleSyms']

        # Loop over dipole terms and find direct products
        dipDictOut = {}
        for s in dipoleDict:
            termsD = terms.copy()
            termsD.append(s)
            # print(termsD)

            # loop over all irreps and check if terms are valid.
            for s2 in self.irreps:
                termsD2 = termsD.copy()
                termsD2.append(s2)
                # print(termsD2)

                rList, rDict = diretProductFromList(PG = self.PG, terms=termsD2)
                # dipDict[s] = rDict.copy()
                dipDictOut[(s,s2)] = dipoleDict[s].copy()
                dipDictOut[(s,s2)]['dipSym'] = s
                dipDictOut[(s,s2)]['targSym'] = s2
                dipDictOut[(s,s2)]['result'] = rList.copy()
                # dipDict[s]['valid'] = [k for k in rDict.keys() if k == 'A1g'] # self.irreps[0]]
                # dipDictOut[(s,s2)]['allowed'] = True if 'A1g' in rList else False
                dipDictOut[(s,s2)]['allowed'] = True if self.irreps[0] in rList else False  # Allowed terms contain totally symmetric rep.
                dipDictOut[(s,s2)]['terms'] = terms
                # dipDict[s].update(dipoleDict[s])

        # Set dataframe and tidy
        pdCont = pd.DataFrame(dipDictOut).drop(['dipSym','targSym']).T
        pdCont.index.rename(['Dipole','Target'], inplace=True)
        pdCont = pdCont.reindex(columns=sorted(pdCont.columns))

        self.continuum = {'dict': dipDictOut.copy(),
                           'PD': pdCont,
                           'allowed': {'PD':pdCont[pdCont.allowed],
                                      'list': pdCont[pdCont.allowed == True ].index.tolist(),   # (Dip, Targ) tuples
                                      'targList': pdCont[pdCont.allowed == True ].index.get_level_values('Target').tolist(),  # Targ only
                                      'dipList': pdCont[pdCont.allowed == True ].index.get_level_values('Dipole').tolist(),}  # Dip only
                          }

        # Display with highlights?
        # For styler see https://www.codespeedy.com/highlight-a-row-in-pandas-data-frame/
        if disp:
            display(pdCont.style.apply(lambda x: ['background-color : yellow']*x.shape[0] if x.allowed else ['background-color : lightgrey']*x.shape[0], axis = 1))   # One liner OK with axis set and mult. by cols.

            
    def scatSym(self, symIon, disp = True):
        """
        Add scattering symmetry column to self.continuum['allowed']['PD'].
        
        scatSym = continuum x ion (per ePS definition)
        
        19/02/24  v1
        
        """
        print(f"Adding scattering symmetries for ion={symIon}"}
        
        scatSymList = []
        
        # Basic version with loop using existing direct product.
        # Should redo as table lookup as directProduct case above?
        for sym in self.continuum['allowed']['PD'].index.get_level_values(level='Target'):
            scatSymList.append(self.directProduct([symIon,sym]))
        
        # Append to main table
        self.continuum['allowed']['PD']['scat'] = scatSymList
        
        if disp:
            display(self.continuum['allowed']['PD'].style.apply(lambda x: ['background-color : yellow']*x.shape[0] if x.allowed else ['background-color : lightgrey']*x.shape[0], axis = 1))   # One liner OK with axis set and mult. by cols.
        

    #*********************** CONVERSION FUNCTIONS

    def toePSproc(self, dimMap = {'C':'Cont','h':'it', 'mu':'muX'}, dataType = 'BLM'):
        """Wrap toePSproc method."""

        self.coeffs[dataType] = toePSproc(self.coeffs, dimMap = dimMap, dataType = dataType, verbose = self.verbose)

    def getIrreps(self):
        """
        Get irrep labels from libmsym repr.

        Quick hack to grab labels from self.ctx.character_table.symmetry_species objects, should be a better way.

        For source, see https://github.com/mcodev31/libmsym/blob/c99470376270db4ec4c925b952fa722e011377d6/bindings/python/libmsym/libmsym.py#L65

        """

        irreps = []

        for item in self.ctx.character_table.symmetry_species:

    #         print(item)
            irreps.append(item.name)

        self.irreps = irreps


    def getSymOps(self):
        """
        Get symmetry operation labels from libmsym repr.

        Quick hack to grab labels from self.ctx.character_table.symmetry_operations objects, should be a better way.

        For source, see https://github.com/mcodev31/libmsym/blob/c99470376270db4ec4c925b952fa722e011377d6/bindings/python/libmsym/libmsym.py#L65

        """

        symOps = []

        for item in self.ctx.character_table.symmetry_operations:

    #         print(item)
            symOps.append(item.__str__().split('(')[1].strip().split(',')[0].split(' ')[0])

        self.symOps = symOps


    def setCharTablePD(self):
        """Generate character table & convert to Pandas DataFrame."""

        charTab = []

        for n,item in enumerate(self.ctx.character_table.symmetry_species):
            charTab.append([item.name, item.dim, *self.ctx.character_table.table[n,:]])
        #     charTab.extend(self.ctx.character_table.table[n,:])
        #     charTab.append([item.name, item.dim])   #.extend(self.ctx.character_table.table[n,:]))

        charTabNP = np.asarray(charTab)

        self.getSymOps()  # Set self.symOps for symmetry operation labels.


        mInd = pd.MultiIndex.from_arrays(charTabNP[:,0:2].T, names=['Character','dim'])
        # df = pd.DataFrame(charTabNP[:,2:], index = mInd, columns=self.ctx.character_table.symmetry_operations[0]._names)
        df = pd.DataFrame(charTabNP[:,2:], index = mInd, columns=self.symOps)

        self.charTablePD = df



    def setCoeffsPD(self, key = 'libmsym', dtype = 'real'):
        """
        Convert raw list output to Pandas DataFrame.

        Parameters
        ----------
        key : str, optional, default = 'libmsym'
            Key for self.coeffs[key]

        dtype : str, optional, default = 'real'
            Key for self.coeffs[key][dtype]

        """

        tabOutNP = self.coeffs[key][dtype]   # For structured array version

        # tabOutNP = np.asarray(self.coeffs['libmsym']['real'])  # Numpy OK, but seems to convert types?
                               # Homogeneous type? Would be OK for values only?
                               # Should be able to set, e.g. table = np.asarray(symObj.coeffTableC) #, dtype='str, int, int, int, int, float')
                               # But currently fails, probably row/col ordering issue? See https://numpy.org/doc/stable/reference/arrays.dtypes.html
#         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character', 'SALC (X)', 'PFIX (h)', 'l', 'm'])

        #********************************** 02/03/22 - updated labels, may break things later!
        # TODO: unify list, use short + long forms?
        # E.g. https://stackoverflow.com/questions/48243818/display-column-name-different-from-dictionary-key-name-in-pandas/48245285
        # Applymap for display? https://pandas.pydata.org/docs/user_guide/style.html#Acting-on-the-Index-and-Column-Headers
        # Or indexing long names: https://stackoverflow.com/questions/49600157/pandas-dataframe-columns-with-long-names
        # Sticky headers: https://pandas.pydata.org/docs/user_guide/style.html#Sticky-Headers
        #
        # UPDATE: set long names to metadata (.attrs['indexes']) and use for display.  See self.displayXlm() for more.
        #
        #**********************************

        #         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character ($\Gamma$)', 'SALC (h)', 'PFIX ($\mu$)', 'l', 'm'])

        # defaultInd = ['C', 'h', 'mu', 'l', 'm']
        defaultInd = self.dims

        # Basic case
        # mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=defaultInd)
        # df = pd.DataFrame(tabOutNP[:,-1].astype(np.float), index = mInd, columns=['b (real)'])

        # For structured array version
        mInd = pd.MultiIndex.from_tuples(tabOutNP[:][self.dims].tolist(), names=defaultInd)  # Set index from named dims
        # df = pd.DataFrame(tabOutNP[:][tabOutNP.dtype.names[-1]], index = mInd, columns=[f'b (dtype)'])  # Set data from final dim in array
        df = pd.DataFrame(tabOutNP[:][tabOutNP.dtype.names[-1]], index = mInd, columns=['b'])

        # mInd
        df.attrs['indexes']= {'shortnames': defaultInd,
                            'longnames': {'C':'Character ($\Gamma$)', 'h':'SALC (h)', 'mu':'PFIX ($\mu$)'},  #, 'l':'l', 'm':'m'},
                            'shortSymbol': {'C':'$\Gamma$', 'mu':'$\mu$'},  #, 'l', 'm'],
                            'libmsym': {'C':'Character', 'h':'SALC', 'mu':'PFIX'},   # 'l', 'm'],
                            'note': "Column name mapping for different sources/conventions. 'shortnames' is default list as used at DF creation."
                            }
        # df.attrs['type'] = 'real'  # REal harmonics from libmsym
        df.attrs['type'] = dtype  # Use passed type

        # TODO: tidy up output & set names & metadata
#         df.name = f"Symmetrize harmonic coeffs {self.PG}"  # Not working?

#         self.coeffDF = df.unstack(level='l').fillna('')  # Set cols by l   - DO THIS ONLY FOR DISPLAY LATER? OR ADD A SWITCH/OPTION?
        # self.coeffDF = df  #.unstack(level='l').fillna('')

        if not 'DF' in self.coeffs.keys():
            self.coeffs['DF'] = {}

        self.coeffs['DF'][dtype] = df
#         .columns.names
        # self.coeffDF.attrs = df.attrs   # Propagate attribs
        # self.coeffs['DF'][dtype].attrs = df.attrs   # Propagate attribs


    def setCoeffsSH(self, absM = True):
        r"""
        Convert symmertrized harmonics to SHtools object, and convert type to complex.

        Parameters
        ----------

        absM : bool, optional, default = True
            Use absM values from input coeff labels?
            May need to force abs(M) for libmsym results? Or double up +/-M terms?

        Notes
        -----

        - Are sign convensions consistent?
           - Looks like +m in libmsym output == symmetric case (same sign for +/-m).
           - And -m is antisym case.
           - Now fixed above (hopefully) by setting -m term as mSign*clmC.coeffs[1,l,m]]
           - Could also be Condon-Shortley phase, have INCLUDED it here (see https://shtools.oca.eu/shtools/public/complex-spherical-harmonics.html#condon-shortley-phase-factor)
        - Current version seems to match Chandra 1987 & Boiko et. al. 2006 for some tested Td cases, more testing required.

        Sources:

        - <div class="csl-bib-body" style="line-height: 1.35; ">
          <div class="csl-entry">Boiko, D.L., Féron, P. and Besnard, P. (2006) ‘Simple method to obtain symmetry harmonics of point groups’, <i>The Journal of Chemical Physics</i>, 125(9), p. 094104. doi:<a href="https://doi.org/10.1063/1.2338040">10.1063/1.2338040</a>.</div>
          <span class="Z3988" title="url_ver=Z39.88-2004&amp;ctx_ver=Z39.88-2004&amp;rfr_id=info%3Asid%2Fzotero.org%3A2&amp;rft_id=info%3Adoi%2F10.1063%2F1.2338040&amp;rft_val_fmt=info%3Aofi%2Ffmt%3Akev%3Amtx%3Ajournal&amp;rft.genre=article&amp;rft.atitle=Simple%20method%20to%20obtain%20symmetry%20harmonics%20of%20point%20groups&amp;rft.jtitle=The%20Journal%20of%20Chemical%20Physics&amp;rft.stitle=J.%20Chem.%20Phys.&amp;rft.volume=125&amp;rft.issue=9&amp;rft.aufirst=D.%20L.&amp;rft.aulast=Boiko&amp;rft.au=D.%20L.%20Boiko&amp;rft.au=P.%20F%C3%A9ron&amp;rft.au=P.%20Besnard&amp;rft.date=2006-09-07&amp;rft.pages=094104&amp;rft.issn=0021-9606"></span>
        </div>



        - <div class="csl-bib-body" style="line-height: 1.35; ">
          <div class="csl-entry">Chandra, N. (1987) ‘Photoelectron spectroscopic studies of polyatomic molecules. I. Theory’, <i>Journal of Physics B: Atomic and Molecular Physics</i>, 20(14), pp. 3405–3415. doi:<a href="https://doi.org/10.1088/0022-3700/20/14/013">10.1088/0022-3700/20/14/013</a>.</div>
          <span class="Z3988" title="url_ver=Z39.88-2004&amp;ctx_ver=Z39.88-2004&amp;rfr_id=info%3Asid%2Fzotero.org%3A2&amp;rft_id=info%3Adoi%2F10.1088%2F0022-3700%2F20%2F14%2F013&amp;rft_val_fmt=info%3Aofi%2Ffmt%3Akev%3Amtx%3Ajournal&amp;rft.genre=article&amp;rft.atitle=Photoelectron%20spectroscopic%20studies%20of%20polyatomic%20molecules.%20I.%20Theory&amp;rft.jtitle=Journal%20of%20Physics%20B%3A%20Atomic%20and%20Molecular%20Physics&amp;rft.volume=20&amp;rft.issue=14&amp;rft.aufirst=N&amp;rft.aulast=Chandra&amp;rft.au=N%20Chandra&amp;rft.date=1987-07-28&amp;rft.pages=3405-3415&amp;rft.spage=3405&amp;rft.epage=3415&amp;rft.issn=0022-3700"></span>
        </div>

        """

        # Set data in
        df = self.coeffs['DF']['real']   # Stacked case
#         df = self.coeffDF.stack(level='l').swaplevel(i='l',j='m')  # TODO: better dim handling here/below?
                                                                # NOTE: swaplevel to force (l,m) order.
#         npTab = np.asarray(self.coeffTable)

#         absM = True  # May need to force abs(M) for libmsym results? Or double up +/-M terms?

        # Get symmetries from Dataframe - assume first dim for this.
        symList = df.index.unique(level=self.dims[0])

        # lmax from Dataframe
        lmax = np.asarray(df.index.unique(level='l')).astype(int).max()

        tabOut = []
        clmSets = {}

        clmInd = 1

        # For each symmetry, get coeffs, push to SHTOOLS format & convert
        for sym in symList:

#             print(sym)

            # Unified version - need to replace empty strings too!
#             npTab = df.loc[sym].droplevel(['SALC (X)', 'PFIX (h)']).reset_index().replace('','NaN').to_numpy().astype('float')  # Original names
            npTab = df.loc[sym].droplevel(['h', 'mu']).reset_index().replace('','NaN').to_numpy().astype('float')  # Short name version - should pull names from attrs or cols however.
            npTab = npTab[~np.isnan(npTab[:,-1])]  # Drop NaNs
#             print(npTab)

            # clm.set_coeffs(npTab)   # NEEDS (values, ls, ms)
            clm = pysh.SHCoeffs.from_zeros(lmax = npTab[:,0].max().astype(int))

            # Store all clms sets (per sym)
#             clmSets[clmInd] = {'sym':sym, 're':clm}
            clmSets[sym] = {'real':clm}

            if absM:
                clm.set_coeffs(npTab[:,2], npTab[:,0].astype(int), np.abs(npTab[:,1].astype(int)))  # OK
            else:
                clm.set_coeffs(npTab[:,2], npTab[:,0].astype(int), npTab[:,1].astype(int))  # OK

        #     clm.coeffs
            clmC = clm.convert(csphase=-1, kind='complex')
#             clmSets[clmInd].update('im':clmC)
#             clmInd += 1
            clmSets[sym].update({'comp':clmC})

            # Tabulate clm entries
            for n in df.loc[sym].index:
            #     print(dfC.loc['A1'].loc[n])

                #*** Single symmetry to numpy array (l,m,coeff)
#                 npTab = df.loc[sym].droplevel(['SALC (X)', 'PFIX (h)']).reset_index().to_numpy()
                #.unstack().to_numpy()
#                 npTab=npTab.astype('float')  # Need to fix type?
            #     npTab
                # Unified version - need to replace empty strings too! NOW SET ABOVE


#                 # clm.set_coeffs(npTab)   # NEEDS (values, ls, ms)
#                 clm = pysh.SHCoeffs.from_zeros(lmax = npTab[:,0].max().astype(int))

#                 if absM:
#                     clm.set_coeffs(npTab[:,2], npTab[:,0].astype(int), np.abs(npTab[:,1].astype(int)))  # OK
#                 else:
#                     clm.set_coeffs(npTab[:,2], npTab[:,0].astype(int), npTab[:,1].astype(int))  # OK

#             #     clm.coeffs
#                 clmC = clm.convert(csphase=-1, kind='complex')


                # Set to output table with updated (l,m) set
                s,h,l,m = n  # Unpack ref values

            #     ind = dfC.loc['A1'].loc[n]

                # Get updated ref value
#                 print(df.loc[n])
        #         print(clm.coeffs[:,int(l),np.abs(int(m))])  # Check inputs
        #         print(clmC.coeffs[:,int(l),np.abs(int(m))]) # Check conversions

                l = int(l)
                if absM:
                    mSign = np.sign(int(m))
                    m = np.abs(int(m))

                else:
                    mSign =1
                    m = int(m)


        #         m = int(m)

        #         tabOut.append([sym, s, h, l, m, clmC.coeffs[:,l,m]])  # This will add 2 values per (l,m)

                # Format +/-m pairs
                try:
                    if m != 0:
                        # tabOut.append([sym, s, h, l, m, clmC.coeffs[0,l,m]])
                        # tabOut.append([sym, s, h, l, -m, mSign*clmC.coeffs[1,l,m]])  # Apply phase fix?
                        tabOut.append((sym, s, h, l, m, clmC.coeffs[0,l,m]))
                        tabOut.append((sym, s, h, l, -m, mSign*clmC.coeffs[1,l,m]))  # Apply phase fix?
                    else:
                        tabOut.append((sym, s, h, l, m, clmC.coeffs[0,l,m]))

                # Skip index errors, can get this is given (l,m) is null/removed for sym group.
                # TODO: should be fixable from PD DataFrame? Due to shared index over all sym groups?
                # NOTE: seemed to be OK, then broke again when playing with index name/styles... bug in index selection somewhere? Should confirm no spurious assignments here.
                except IndexError:
                    pass

        # Set new output tables
#         tabOutNP = np.asarray(tabOut)  # Numpy OK, but seems to convert types?
#                                        # Homogeneous type? Would be OK for values only?
#
# #         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character', 'SALC (X)', 'PFIX (h)', 'l', 'm'])
# #         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character ($\Gamma$)', 'SALC (h)', 'PFIX ($\mu$)', 'l', 'm'])
#         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=self.coeffs['DF']['real'].attrs['indexes']['shortnames'])
#         dfC = pd.DataFrame(tabOutNP[:,-1].astype(complex), index = mInd, columns=['b (complex)'])

        # # For structured array version
        # dimTypes = list(zip(self.dims,self.dtypes))  # Dims
        # dimTypes.extend([('b',complex)])   # Data
        # defaultInd = self.dims
        # tabOutNP = np.asarray(tabOut, dtype=dimTypes)
        # mInd = pd.MultiIndex.from_tuples(tabOutNP[:][self.dims].tolist(), names=defaultInd)  # Set index from named dims
        # dfC = pd.DataFrame(tabOutNP[:][tabOutNP.dtype.names[-1]], index = mInd, columns=['b (complex)'])  # Set data from final dim in array


        # Store outputs
        # self.coeffTableC = tabOut    # Use raw results here?
        # self.coeffs['libmsym']['comp'] = tabOut
        # self.coeffs['DF']['comp'] = dfC
        # self.coeffs['DF']['comp'].attrs = self.coeffs['DF']['real'].attrs  # Use inital attribs
        # self.coeffs['DF']['comp'].attrs['type'] = 'comp'
#         self.clm = {'re':clm, 'im':clmC, 'note':'SHtools clm coefficient objects, re (real) and im (complex) harmonic expansions.'}
        # self.clm = clmSets
        # self.clm.update({'note':'SHtools clm coefficient objects, real and comp (complex) harmonic expansions.'})
        self.coeffs['SH'] = clmSets
        self.coeffs['SH'].update({'note':'SHtools clm coefficient objects, real and comp (complex) harmonic expansions.'})

        # With unified methods
        # Convert to NP structured array (keeps heterogeneous dtypes)
        dimTypes = list(zip(self.dims,self.dtypes))  # Dims
        dimTypes.extend([('b',complex)])   # Data
        self.coeffs['libmsym'] = {'comp':np.asarray(tabOut, dtype=dimTypes)}
        self.setCoeffsPD(dtype = 'comp')


    def setCoeffsXR(self, stack = True):
        """
        Convert coeffs from PD DataFrame to Xarray Dataset.

        Parameters
        ----------

        stack : dict, bool, optional, default = True
            If True, try and use default stacking, {'inds':self.dims[1:3], 'LM':self.dims[3:]}.
            If False, don't stack.
            If dict, try and stack with specified dict mapping.

        """

        if xrFlag:
            coeffsXR = xr.Dataset.from_dataframe(self.coeffs['DF']['real']).rename({'b':'b (real)'})
            coeffsXR.update(xr.Dataset.from_dataframe(self.coeffs['DF']['comp']).rename({'b':'b (comp)'}))

            self.coeffs['XR'] = coeffsXR

            # Set additional metadata for portability
            # TODO: check ePSproc conventions here.
            self.coeffs['XR'].attrs = {'dataType':'symHarm',
                                        'name':'Symmetrized harmonics',
                                        'PG': self.PG,
                                        'lmax': self.lmax}

            # May need to fix coord types from str to ints
            # TODO: better solution here? Should fix in Pandas tables?
            for k in self.coeffs['XR'].coords.keys():
                try:
                    self.coeffs['XR'].coords[k].data = self.coeffs['XR'].coords[k].data.astype(int)
                except ValueError:
                    pass

            if stack:
                if isinstance(stack,dict):
                    self.coeffs['XR'] = self.coeffs['XR'].stack(stack)
                else:
                    self.coeffs['XR'] = self.coeffs['XR'].stack({'inds':self.dims[1:3], 'LM':self.dims[3:]})

        else:
            print("Xarray required to run self.setCoeffsXR.")




    #*********************** DISPLAY FUNCTIONS

    def printCharacterTable(self, returnPD = False):
        """
        Print character table with species & degen using Pandas

        Parameters
        ----------

        returnPD : bool, optional, default = False
            Return PD object instead of display() if True.

        Returns
        -------

        Empty unless returnPD = True set.

        """

        # Basic version
    #         for n,item in enumerate(self.ctx.character_table.symmetry_species):
    #             print(item.name, item.dim, self.ctx.character_table.table[n,:])

        # PD version
        if not hasattr(self, 'charTablePD'):
            self.setCharTablePD()

        if returnPD:
            return self.charTablePD

        # TODO - add isnotebook and related checks here
        display(self.charTablePD)



    # Display table
    def displayXlm(self, names = 'longnames', YlmType = 'real', setCols = 'l',
                    dropLevels=[], returnPD = False, sticky = False,
                    symFilter = False, symFilterChannel = 'Target'):
        """
        Print table of values from Pandas Dataframe self.coeffs['DF']['real'], with specified labels (from self.coeffs['DF']['real'].attrs['indexes']).

        Parameters
        ----------

        names : str, optional, default = 'longnames'
            Labels to use in printed table, from self.coeffs['DF']['real'].attrs['indexes']

        YlmType : str, optional, default = 'real'
            - 'real' show real harmonic coeffs, from self.coeffs['DF']['real']
            - 'comp' show complex harmonic coeffs, from self.coeffs['DF']['comp']

        setCols : str, optional, default = 'l'
            Set which label to use for display.
            Set via self.coeffsDF.unstack(level=setCols).
            Note this level must be in the DataFrame index, and there is currently no error checking.

        dropLevels : str or list of strings, default = []
            Drop levels specified from displayed table.

        returnPD : bool, optional, default = False
            Return PD object instead of display() if True.

        sticky : bool, optional, default = False
            Apply "sticky" index styler to displayed table.
            (If supported by Pandas version.)

        symFilter : bool, optional, default = False
            If True, apply filter from symFilterChannel.

        symFilterChannel : str, optional, default = None
            Try and filter output based on allowed symmetries.
            This requires self.continuum to be set, and to match a column name for filtering.
            E.g. 'Target'  will filter the Xlm table from self.continuum['allowed']['PD']['Target'] symmetries.

            TODO: allow use of short names here, currently only filters on output names.


        Returns
        -------

        Empty unless returnPD = True set.

        """

#         display(self.coeffDF.rename(index=self.coeffDF.attrs['indexes']['longnames']))  # Works for col or row values, not headers (==names)

#         renamed = symObj.coeffDF.index.set_names(symObj.coeffDF.attrs['indexes']['longnames'])  # Not working in Pandas v1.1 .... but should be OK later: https://pandas.pydata.org/docs/reference/api/pandas.Index.set_names.html

        # Manual version works OK - in Pandas v1.1 rename needs full list
#         newNames = [v for k,v in self.coeffDF.attrs['indexes'][names].items() if k in self.coeffDF.index.names]  # Map dictionary

        # Current strucutre - may want to set to dicts?
        # Note hard-coded unstack too - should add as an option
        if YlmType == 'real':
            inputData = self.coeffs['DF']['real'].copy().unstack(level=setCols).fillna('')
            inputData.attrs = self.coeffs['DF']['real'].attrs
        elif YlmType == 'comp':
            inputData = self.coeffs['DF']['comp'].copy().astype(str).unstack(level=setCols).fillna('')  # Note cast to string to avoid `TypeError: No matching signature found` for complex data on unstack (might be PD version bug, tested in v1.0.1)
            inputData.attrs = self.coeffs['DF']['comp'].attrs
        else:
            print(f"Didn't recognise Ylm type {YlmType}")

        # Drop levels if specified (will return as is if dropLevels empty)
        inputData = inputData.droplevel(dropLevels)

        # With variable len remap
        newNames = []
        for k in inputData.index.names:

            if k in inputData.attrs['indexes'][names].keys():
                newNames.append(inputData.attrs['indexes'][names][k])
            else:
                newNames.append(k)

#         test = self.coeffDF.copy()   # Display as is
        test = inputData   # With unstack and fillna
        test.index.rename(names = newNames, inplace=True)

        if symFilter:
            if hasattr(self, 'continuum'):
                subsetMask = test.index.get_level_values(test.index.names[0]).isin(self.continuum['allowed']['PD'].index.get_level_values(symFilterChannel))
                # subsetMask = test.index.get_level_values(newNames[0]).isin(self.continuum['allowed']['PD'].index.get_level_values(symFilter))
                test = test[subsetMask]
            else:
                print(f"*** Cannot filter on symmetry without self.continuum, run self.directProductContinuum() to generate.")


        if returnPD:
            return test

        if sticky:
            # Try displaying with sticky index, but may fail in older PD versions (OK in v1.4, fails in v1.1)
            try:
                display(test.style.set_sticky(axis="index"))
            except AttributeError:
                display(test)

        else:
            display(test)


    def plotXlm(self, pType = 'real', gridlmax = 20, syms = None, **kwargs):
        """
        Quick plot of Xlm by symmetry using SHtools grid.plot().

        Parameters
        ----------
        pType : str, optional, default = 'real'
            Plot type as key (for self.coeffs['SH'][pType]).
            Default cases are 'real' or 'comp' types (complex harmonics)

        gridlmax : int or None, optional, default = 20
            Used by SHtools clm.expand() routine to define gridding.
            Use SHtools defaults if set to None (== lmax of distribution)

       syms : str or list, default = None
           Symmetry groups self.coeffs['SH'][pType][sym] to plot.
           Defaults to all syms, as defined by self.coeffs['SH'][pType].keys()


        **kwargs : optional args passed to SHtools grid.plot(**kwargs)

        Notes
        -----

        - Gridding is defined automatically by clm.expand() routine, pass lmax=int as a proxy to override defaults


        TODO: For more control with subselection etc. needs bridging with other plotting tools (BLMplot etc.).

        """
        # Set syms
        if syms is None:
            syms = self.coeffs['SH'].keys()

        if type(syms) is not list:
            syms = list(syms)

        # Plot
        for key in syms:
            if key != 'note':
#                 print(key)
                grid = self.coeffs['SH'][key][pType].expand(lmax = gridlmax)

                # 'title' fails for some versions of SHtools/MatPlotlib?
                # BUT this code shows blank plot for try part, even on fail, and with show=False set.
                # try:
                #     grid.plot(title=f'{self.PG} ({key}), l_max={self.lmax}', show=False, **kwargs)
                # except AttributeError:
                #     grid.plot(**kwargs)

                # More general fix, parse p and set titles (this is (fig, [axes]) object)
                p = grid.plot(**kwargs)
                # p[1].set_title('test');
                prefixAxisTitles(p, prefix = f'{self.PG} ({key}), l_max={self.lmax}')  # Use separate function to parse p and set title(s)
                                                                                        # This works also for re/im case, with nested plots.

#                 fig, ax = grid.plot(show=False)
        #         fig.set_label(key)
        #         ax.title = key
