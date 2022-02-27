"""
Class for determining and handling symmetrized harmonics.

24/02/22    v1 in development.

- For early dev work see http://localhost:8888/lab/workspaces/symm/tree/python/chem/tools/symmetrized_harmonics_libmsym_tests_160122.ipynb
- For class dev see http://localhost:8888/lab/workspaces/symm/tree/python/chem/tools/symmetrized_harmonics_PEMtk-dev_240222.ipynb

"""

import libmsym as msym

import pyshtools as pysh

import numpy as np
import pandas as pd

import argparse, random, sys

class symHarm():
    """
    Class to set and manipulate symmetrized harmonics.

    Main symmetry routine uses [libmsym](https://github.com/mcodev31/libmsym), code adapted from https://github.com/mcodev31/libmsym/issues/21

    """

    def __init__(self, PG = 'Cs', lmax = 4):

        # Set params
        self.PG = PG
        self.lmax = lmax

        # Set single element as expansion centre
        self.elements = [msym.Element(name = "C", coordinates = [0.0, 0.0,0.0])]

        # Compute harmonic coeffs
        self.set_lm_basis()
        self.calcSymHarmonics()


    def set_lm_basis(self):
        """Set basis functions for each (element, l, m)."""

        self.basis_functions = []

        for l in range(0, self.lmax):
            for m in range(-l,l+1):
        #         basis_functions.extend([set_basis(e, l=l, m=m, n=n, name=f"{l}s{m}") for e in elements])
                self.basis_functions.extend([self.set_basis(e, l=l, m=m, n=self.lmax, name=f"{l},{m}") for e in self.elements])


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


    def printCharacterTable(self):
        """Print character table with species & degen"""

        # Basic version
#         for n,item in enumerate(self.ctx.character_table.symmetry_species):
#             print(item.name, item.dim, self.ctx.character_table.table[n,:])

        # PD version
        if not hasattr(self, 'charTablePD'):
            self.setCharTablePD()

        # TODO - add isnotebook and related checks here
        display(self.charTablePD)


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
        mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character', 'SALC (X)', 'PFIX (h)', 'l', 'm'])
        df = pd.DataFrame(tabOutNP[:,-1], index = mInd, columns=['Value'])
        # mInd

        # TODO: tidy up output & set names & metadata
#         df.name = f"Symmetrize harmonic coeffs {self.PG}"  # Not working?

        self.coeffDF = df.unstack(level='l').fillna('')  # Set cols by l
#         .columns.names



    def setCoeffsSH(self, absM = True):
        """
        Convert symmertrized harmonics to SHtools object, and convert type to complex.

        Parameters
        ----------

        absM : bool, optional, default = True
            Use absM values from input coeff labels?
            May need to force abs(M) for libmsym results? Or double up +/-M terms?

        """

        # Set data in
        df = self.coeffDF.stack(level='l').swaplevel(i='l',j='m')  # TODO: better dim handling here/below?
                                                                # NOTE: swaplevel to force (l,m) order.
#         npTab = np.asarray(self.coeffTable)

#         absM = True  # May need to force abs(M) for libmsym results? Or double up +/-M terms?

        # Get symmetries from Dataframe
        symList = df.index.unique(level='Character')

        # lmax from Dataframe
        lmax = np.asarray(df.index.unique(level='l')).astype(int).max()

        tabOut = []
        clmSets = {}

        clmInd = 1

        # For each symmetry, get coeffs, push to SHTOOLS format & convert
        for sym in symList:

#             print(sym)

            # Unified version - need to replace empty strings too!
            npTab = df.loc[sym].droplevel(['SALC (X)', 'PFIX (h)']).reset_index().replace('','NaN').to_numpy().astype('float')
            npTab = npTab[~np.isnan(npTab[:,-1])]  # Drop NaNs
#             print(npTab)

            # clm.set_coeffs(npTab)   # NEEDS (values, ls, ms)
            clm = pysh.SHCoeffs.from_zeros(lmax = npTab[:,0].max().astype(int))

            # Store all clms sets (per sym)
#             clmSets[clmInd] = {'sym':sym, 're':clm}
            clmSets[sym] = {'re':clm}

            if absM:
                clm.set_coeffs(npTab[:,2], npTab[:,0].astype(int), np.abs(npTab[:,1].astype(int)))  # OK
            else:
                clm.set_coeffs(npTab[:,2], npTab[:,0].astype(int), npTab[:,1].astype(int))  # OK

        #     clm.coeffs
            clmC = clm.convert(csphase=-1, kind='complex')
#             clmSets[clmInd].update('im':clmC)
#             clmInd += 1
            clmSets[sym].update({'im':clmC})

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
        #         print(df.loc[n])
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
#                 try:
                if m != 0:
                    tabOut.append([sym, s, h, l, m, clmC.coeffs[0,l,m]])
                    tabOut.append([sym, s, h, l, -m, mSign*clmC.coeffs[1,l,m]])  # Apply phase fix?
                else:
                    tabOut.append([sym, s, h, l, m, clmC.coeffs[0,l,m]])

                # Skip index errors, can get this is given (l,m) is null/removed for sym group.
                # TODO: should be fixable from PD DataFrame? Due to shared index over all sym groups?
#                 except IndexError:
#                     pass

        # Set new output tables
        tabOutNP = np.asarray(tabOut)  # Numpy OK, but seems to convert types?
                                       # Homogeneous type? Would be OK for values only?


        mInd = pd.MultiIndex.from_arrays(tabOutNP[:,:-1].T, names=['Character', 'SALC (X)', 'PFIX (h)', 'l', 'm'])
        dfC = pd.DataFrame(tabOutNP[:,-1], index = mInd, columns=['clmC'])


        # Store outputs
        self.coeffTableC = tabOut    # Use raw results here?
        self.coeffDFC = dfC
#         self.clm = {'re':clm, 'im':clmC, 'note':'SHtools clm coefficient objects, re (real) and im (complex) harmonic expansions.'}
        self.clm = clmSets
        self.clm.update({'note':'SHtools clm coefficient objects, re (real) and im (complex) harmonic expansions.'})
