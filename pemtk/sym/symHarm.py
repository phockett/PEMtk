"""
Class for determining and handling symmetrized harmonics.

24/02/22    v1 in development.

- For early dev work see http://localhost:8888/lab/workspaces/symm/tree/python/chem/tools/symmetrized_harmonics_libmsym_tests_160122.ipynb
- For class dev see http://localhost:8888/lab/workspaces/symm/tree/python/chem/tools/symmetrized_harmonics_PEMtk-dev_240222.ipynb

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
    - Spherical harmonic expansions and conversions (real, imaginary, normalization etc.) and basic ploting (2D maps) are handled by [pySHtools](https://shtools.oca.eu).
    - TODO: Interface to PEMtk/ePSproc for other plotters and handling routines.


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


    def __init__(self, PG = 'Cs', lmax = 4):
        """Init class object"""

        # Set params
        self.PG = PG
        self.lmax = lmax

        # Set single element as expansion centre
        self.elements = [msym.Element(name = "C", coordinates = [0.0, 0.0,0.0])]

        # Compute harmonic coeffs
        self.set_lm_basis()
        self.calcSymHarmonics()
        self.setCoeffsSH()


    #*********************** CALCULATION FUNCTIONS

    def set_lm_basis(self):
        """Set basis functions for each (element, l, m)."""

        self.basis_functions = []

        for l in range(0, self.lmax+1):
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
                            tabOut.append([ctx.character_table.symmetry_species[ss.symmetry_species].name, salcix, pfix, int(l), int(m), v])

        # Log outputs - may want to use dicts?
        self.ctx = ctx
        self.coeffTable = tabOut

        # Set PD table form too
        self.setCoeffsPD()




    #*********************** CONVERSION FUNCTIONS

    def setCharTablePD(self):
        """Generate character table & convert to Pandas DataFrame."""

        charTab = []

        for n,item in enumerate(self.ctx.character_table.symmetry_species):
            charTab.append([item.name, item.dim, *self.ctx.character_table.table[n,:]])
        #     charTab.extend(self.ctx.character_table.table[n,:])
        #     charTab.append([item.name, item.dim])   #.extend(self.ctx.character_table.table[n,:]))

        charTabNP = np.asarray(charTab)


        mInd = pd.MultiIndex.from_arrays(charTabNP[:,0:2].T, names=['Character','dim'])
        df = pd.DataFrame(charTabNP[:,2:], index = mInd, columns=symObj.ctx.character_table.symmetry_operations[0]._names)

        self.charTablePD = df



    def setCoeffsPD(self):
        """Convert raw list output to Pandas DataFrame,"""

        tabOutNP = np.asarray(self.coeffTable)  # Numpy OK, but seems to convert types?
                               # Homogeneous type? Would be OK for values only?
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

        defaultInd = ['C', 'h', 'mu', 'l', 'm']
        mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=defaultInd)
        df = pd.DataFrame(tabOutNP[:,-1], index = mInd, columns=['b (real)'])
        # mInd
        df.attrs['indexes']= {'shortnames': defaultInd,
                            'longnames': {'C':'Character ($\Gamma$)', 'h':'SALC (h)', 'mu':'PFIX ($\mu$)'},  #, 'l':'l', 'm':'m'},
                            'shortSymbol': {'C':'$\Gamma$', 'mu':'$\mu$'},  #, 'l', 'm'],
                            'libmsym': {'C':'Character', 'h':'SALC', 'mu':'PFIX'},   # 'l', 'm'],
                            'note': "Column name mapping for different sources/conventions. 'shortnames' is default list as used at DF creation."
                            }
        df.attrs['type'] = 'real'  # REal harmonics from libmsym

        # TODO: tidy up output & set names & metadata
#         df.name = f"Symmetrize harmonic coeffs {self.PG}"  # Not working?

#         self.coeffDF = df.unstack(level='l').fillna('')  # Set cols by l   - DO THIS ONLY FOR DISPLAY LATER? OR ADD A SWITCH/OPTION?
        self.coeffDF = df  #.unstack(level='l').fillna('')
#         .columns.names
        self.coeffDF.attrs = df.attrs   # Propagate attribs


    def setCoeffsSH(self, absM = True):
        """
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
        df = self.coeffDF   # Stacked case
#         df = self.coeffDF.stack(level='l').swaplevel(i='l',j='m')  # TODO: better dim handling here/below?
                                                                # NOTE: swaplevel to force (l,m) order.
#         npTab = np.asarray(self.coeffTable)

#         absM = True  # May need to force abs(M) for libmsym results? Or double up +/-M terms?

        # Get symmetries from Dataframe
        symList = df.index.unique(level='C')

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
                        tabOut.append([sym, s, h, l, m, clmC.coeffs[0,l,m]])
                        tabOut.append([sym, s, h, l, -m, mSign*clmC.coeffs[1,l,m]])  # Apply phase fix?
                    else:
                        tabOut.append([sym, s, h, l, m, clmC.coeffs[0,l,m]])

                # Skip index errors, can get this is given (l,m) is null/removed for sym group.
                # TODO: should be fixable from PD DataFrame? Due to shared index over all sym groups?
                # NOTE: seemed to be OK, then broke again when playing with index name/styles... bug in index selection somewhere? Should confirm no spurious assignments here.
                except IndexError:
                    pass

        # Set new output tables
        tabOutNP = np.asarray(tabOut)  # Numpy OK, but seems to convert types?
                                       # Homogeneous type? Would be OK for values only?


#         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character', 'SALC (X)', 'PFIX (h)', 'l', 'm'])
#         mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character ($\Gamma$)', 'SALC (h)', 'PFIX ($\mu$)', 'l', 'm'])
        mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=self.coeffDF.attrs['indexes']['shortnames'])
        dfC = pd.DataFrame(tabOutNP[:,-1], index = mInd, columns=['b (complex)'])


        # Store outputs
        self.coeffTableC = tabOut    # Use raw results here?
        self.coeffDFC = dfC
        self.coeffDFC.attrs = self.coeffDF.attrs  # Use inital attribs
        self.coeffDFC.attrs['type'] = 'comp'
#         self.clm = {'re':clm, 'im':clmC, 'note':'SHtools clm coefficient objects, re (real) and im (complex) harmonic expansions.'}
        self.clm = clmSets
        self.clm.update({'note':'SHtools clm coefficient objects, real and comp (complex) harmonic expansions.'})


    def setCoeffsXR(self):
        """Convert coeffs from PD DataFrame to Xarray Dataset"""

        if xrFlag:
            coeffsXR = xr.Dataset.from_dataframe(self.coeffDF)
            coeffsXR.update(xr.Dataset.from_dataframe(self.coeffDFC))

            self.coeffsXR = coeffsXR

        else:
            print("Xarray required to run self.setCoeffsXR.")




    #*********************** DISPLAY FUNCTIONS

    def printCharacterTable(self):
        """Print character table with species & degen using Pandas"""

        # Basic version
    #         for n,item in enumerate(self.ctx.character_table.symmetry_species):
    #             print(item.name, item.dim, self.ctx.character_table.table[n,:])

        # PD version
        if not hasattr(self, 'charTablePD'):
            self.setCharTablePD()

        # TODO - add isnotebook and related checks here
        display(self.charTablePD)



    # Display table
    def displayXlm(self, names = 'longnames', YlmType = 'real'):
        """
        Print table of values from self.coeffDF, with specified labels (from self.coeffDF.attrs['indexes']).

        Parameters
        ----------

        names : str, optional, default = 'longnames'
            Labels to use in printed table, from self.coeffDF.attrs['indexes']

        YlmType : str, optional, default = 'real'
            - 'real' show real harmonic coeffs, from self.coeffDF
            - 'comp' show complex harmonic coeffs, from self.coeffDFC

        """

#         display(self.coeffDF.rename(index=self.coeffDF.attrs['indexes']['longnames']))  # Works for col or row values, not headers (==names)

#         renamed = symObj.coeffDF.index.set_names(symObj.coeffDF.attrs['indexes']['longnames'])  # Not working in Pandas v1.1 .... but should be OK later: https://pandas.pydata.org/docs/reference/api/pandas.Index.set_names.html

        # Manual version works OK - in Pandas v1.1 rename needs full list
#         newNames = [v for k,v in self.coeffDF.attrs['indexes'][names].items() if k in self.coeffDF.index.names]  # Map dictionary

        # Current strucutre - may want to set to dicts?
        # Note hard-coded unstack too - should add as an option
        if YlmType is 'real':
            inputData = self.coeffDF.copy().unstack(level='l').fillna('')
            inputData.attrs = self.coeffDF.attrs
        elif YlmType is 'comp':
            inputData = self.coeffDFC.copy().unstack(level='l').fillna('')
            inputData.attrs = self.coeffDFC.attrs
        else:
            print(f"Didn't recognise Ylm type {YlmType}")

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

        # Try displaying with sticky index, but may fail in older PD versions (OK in v1.4, fails in v1.1)
        try:
            display(test.style.set_sticky(axis="index"))
        except AttributeError:
            display(test)


    def plotXlm(self, pType = 'real', gridlmax = 20, syms = None, **kwargs):
        """
        Quick plot of Xlm by symmetry using SHtools grid.plot().

        Parameters
        ----------
        pType : str, optional, default = 'real'
            Plot type as key (for self.clm[pType]).
            Default cases are 'real' or 'comp' types (complex harmonics)

        gridlmax : int or None, optional, default = 20
            Used by SHtools clm.expand() routine to define gridding.
            Use SHtools defaults if set to None (== lmax of distribution)

       syms : str or list, default = None
           Symmetry groups self.clm[pType][sym] to plot.
           Defaults to all syms, as defined by self.clm[pType].keys()


        **kwargs : optional args passed to SHtools grid.plot(**kwargs)

        Notes
        -----

        - Gridding is defined automatically by clm.expand() routine, pass lmax=int as a proxy to override defaults


        TODO: For more control with subselection etc. needs bridging with other plotting tools (BLMplot etc.).

        """
        # Set syms
        if syms is None:
            syms = self.clm.keys()

        if type(syms) is not list:
            syms = list(syms)

        # Plot
        for key in syms:
            if key != 'note':
#                 print(key)
                grid = self.clm[key][pType].expand(lmax = gridlmax)
                grid.plot(title=f'{self.PG} ({key}), l_max={self.lmax}', **kwargs)

#                 fig, ax = grid.plot(show=False)
        #         fig.set_label(key)
        #         ax.title = key
