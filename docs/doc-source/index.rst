.. PEMtk
   Adapted from ePSproc documentation master file, created by
   sphinx-quickstart on Tue Aug  6 15:55:35 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   06/02/20 - Updated with nbsphinx support.

Welcome to PEMtk - the Photoelectron Metrology Toolkit
======================================================

(aka Quantum Metrology with Photoelectrons platform data & analysis
layer: the unifiying layer (glue) for the platform)

Currently, the platform is under development, based on methods developed over the last 10+ years of reasearch in this area. On the theory side, work is based around ePolyScat, and a set of python packages have already been developed; on the experimental side, the plan is to update existing Matlab codes for Velocity Map Imaging (VMI) experiments and analysis routines (and rewrite/unify in python). The real foundation, and glue, for the platform will be the Photoelectron metrology toolkit (PEMtk), which will provide the unifying data platform, and analysis routines. In the future, it is hoped that this platform will be extended to other theoretical and experimental methods, but continue to provide a useful, unifying, platform.

See the `intro for more details <about.html>`__.

.. toctree::
   :maxdepth: 2
   :caption: Intro:

   about

.. toctree::
   :maxdepth: 2
   :caption: Symmetrized harmonics:

   sym/pemtk_symHarm_demo_160322_tidy
   sym/pemtk_symHarm_epsproc-interface_demo_240322


.. toctree::
   :maxdepth: 2
   :caption: Fitting/extracting matrix elements:

   fitting/PEMtk_fitting_basic_demo_030621-full_010922
   .. fitting/PEMtk_fitting_demo_multi-fit_tests_130621_r1000
   fitting/PEMtk_fitting_demo_multi-fit_tests_130621-para_010922
   fitting/PEMtk_analysis_demo_150621-tidy
   fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy
   fitting/PEMtk_fitting_basis-set_demo_050621-full
   fitting/PEMtk_fitting_backends_demo_010922


.. toctree::
   :maxdepth: 4
   :caption: Function ref:

   modules/pemtk


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
