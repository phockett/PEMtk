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
