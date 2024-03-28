PEMtk: the Photoelectron Metrology Toolkit
==========================================

|zenodo|

Toolkit for photoelectron metrology, data & analysis layer.

Status: as of late 2023, full modelling and fitting of LF, AF and MF observables is implemented, and released as v0.0.1 of the package. Note this is still beta code however, and in active development.

For details see:

1. `Ongoing PEMtk docs on ReadTheDocs <https://pemtk.readthedocs.io/en/latest/index.html>`__ - includes getting started guide and basic demos.
2. `Quantum Metrology with Photoelectrons Vol. 3 <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/intro.html>`__ - open source book covering the theory and numerical implementation in detail, including Jupyter notebooks and full numerical examples of matrix element retrieval for three demonstration cases (N2, OCS, C2H4), using `v0.0.1-dev release <https://github.com/phockett/PEMtk/releases/tag/v0.0.1-dev-QM3-310723>`__ of the codebase.
3. `Topical Review: Extracting Molecular Frame Photoionization Dynamics from Experimental Data (2022) <https://www.authorea.com/users/71114/articles/447808-extracting-molecular-frame-photoionization-dynamics-from-experimental-data>`__ - manuscript which covers photoionization theory and fitting methodologies, the latter as somewhat of a manual for PEMtk, including a full numerical case study.

The figures below illustrate results for matrix element retrieval from N2, based on simulated aligned-frame data, illustrating density matrix and MFPAD results and phase sensitivity. For more details, see the docs above, or [the source notebooks](https://pemtk.readthedocs.io/en/latest/topical_review_case_study/matrix_element_extraction_MFrecon_PEMtk_180722-dist.html).

.. figure:: https://raw.githubusercontent.com/phockett/PEMtk/0a40bf2b38cff8187b2265094b4d7d0e8c8ee17e/docs/doc-source/figs/MFPADs_N2_recon_demo_2023.png
  :alt: MFPAD reconstruction demo

  MFPAD reconstruction demo for N2, see `Quantum Metrology with Photoelectrons Vol. 3 <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/part2/case-study-N2_290723.html#plot-mf-pads>`__ for details.


.. figure:: https://raw.githubusercontent.com/phockett/PEMtk/0a40bf2b38cff8187b2265094b4d7d0e8c8ee17e/docs/doc-source/figs/denMat_N2_recon_demo_2023.png
  :alt: Density matrix reconstruction demo

  Density matrix reconstruction demo for N2, see `Quantum Metrology with Photoelectrons Vol. 3 <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/part2/case-study-N2_290723.html#density-matrices>`__ for details.



Currently implemented (v0.0.1)
------------------------------

- `Symmetrized harmonics generation functions (March 2022) <https://pemtk.readthedocs.io/en/latest/sym/pemtk_symHarm_demo_160322_tidy.html>`__. (Requires `libmsym <https://github.com/mcodev31/libmsym>`__.)
- `Fitting analysis functions and plotters (Dec. 2021) <https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy.html>`__. (Requires `lmfit <https://lmfit.github.io/lmfit-py/intro.html>`__.)
- Basic data model (from `ePSproc base class <https://epsproc.readthedocs.io/en/latest/demos/ePSproc_class_demo_161020.html>`__ currently) and `fitting class (20/05/21) + docs <https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html>`__. (Requires `ePSproc <https://epsproc.readthedocs.io>`__.)
- Code framework and notes (August 2020).

A full demonstration can be found in `Quantum Metrology with Photoelectrons Vol. 3 (2023) <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/intro.html>`__, which made use of the `v0.0.1-dev release of the codebase (Aug 2023) <https://github.com/phockett/PEMtk/releases/tag/v0.0.1-dev-QM3-310723>`__.


Installation
------------

Currently (dev version) only from source + local pip install

Using pip

.. code-block::

  git clone https://github.com/phockett/PEMtk.git
  pip install -e PEMtk


Or from setup.py (from clone dir)

``python setup.py install``



Notes

* The dev version doesn't force any dependencies upon installation, but most functionality depends on `ePSproc` and required packages, which can be installed via pip or manually, see `ePSproc installation notes <https://epsproc.readthedocs.io/en/latest/about.html#installation-python>`__. Some optional packages may require further installation efforts, see the `scripts included in the Docker builds for recipes <https://github.com/phockett/open-photoionization-docker-stacks/tree/main/epsproc-pemtk>`__.
* The repo can be passed directly to pip, e.g. ``pip install git+https://github.com/phockett/PEMtk.git``, see `notes in the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#git>`_.
* Note that ``pip -e`` is for 'editable', and requires the source dir to remain, but the installation is also editable, `see notes here <https://stackoverflow.com/questions/41535915/python-pip-install-from-local-dir>`_. Alternatively, use ``pip install <path_to_local_pkg>``.


TODO: tidy-up build/install chain.


Docker
------

Docker builds are available as follows:

- `Full image build via Dockerhub for the QM3 book (Aug 2023) <https://hub.docker.com/r/epsproc/quantum-met-vol3>`__. (Also available `via Zenodo DOI:10.5281/zenodo.8286020 <https://doi.org/10.5281/zenodo.8286020>`__.)
  - To use this image, simply run `docker pull epsproc/quantum-met-vol3` to pull a copy, then `docker run epsproc/quantum-met-vol3` to run with default settings (which uses port 8888 for JupyterLab). The Jupyter Lab interface will be available at http://localhost:8888, with default password `qm3`.
  - For further details see the `QM3 Docker builds notes <https://github.com/phockett/Quantum-Metrology-with-Photoelectrons-Vol3#docker-builds>`__.
- Source Dockerfiles for ePSproc + PEMtk builds from `the Open Photoionization Docker Stacks <https://github.com/phockett/open-photoionization-docker-stacks/tree/main/epsproc-pemtk>`__.


Roadmap
-------

- Further fitting methodology & code developments
  - Faster fitting, GPU fitting etc.
  - Implementation of matrix techniques (see Gregory, M. et al. (2021) ‘Towards molecular frame photoelectron angular distributions in polyatomic molecules from lab frame coherent rotational wavepacket evolution’, Journal of Physics B: Atomic, Molecular and Optical Physics, 54(14), p. 145601. Available at: https://doi.org/10.1088/1361-6455/ac135f.)
- Integrated data handling class(es).
- Image processing (basic inversions, tomography, FT methods etc.; some aspects are already implemented in `TMO-dev package <https://github.com/phockett/tmo-dev>`__).
- Image simulation.
- More versatile plotting routines.
- Interfaces for various experimental platforms.


Quantum Metrology with Photoelectrons Platform
----------------------------------------------

PEMtk, aka the Quantum Metrology with Photoelectrons platform data & analysis layer, is the unifying layer (glue) for the platform.

Currently, the platform is under development, based on methods developed over the last 10+ years of research in this area. On the theory side, work is based around ePolyScat, and a set of python packages have already been developed (`ePSproc <https://epsproc.readthedocs.io>`__); on the experimental side, the plan is to update existing Matlab codes for Velocity Map Imaging (VMI) experiments and analysis routines (and rewrite/unify in python). Some early work in python - specifically for FEL data - can be found in the `TMO-dev package <https://github.com/phockett/tmo-dev>`__. The real foundation, and glue, for the platform will be the Photoelectron metrology toolkit (PEMtk), which will provide the unifying data platform, and analysis routines. In the future, it is hoped that this platform will be extended to other theoretical and experimental methods, but continue to provide a useful, unifying, platform. See also `ePSdata for general aims & motivation <https://phockett.github.io/ePSdata/about.html#Motivation>`__, and a growing collection of computational results for various cases.

.. Local fig: .. figure:: ./docs/doc-source/figs/QM_unified_schema_wrapped_280820.gv.png
   Use GH version via full URL instead for consistency on RTD.

.. figure:: https://raw.githubusercontent.com/phockett/PEMtk/4eec9217203bfd1aee13bd8b64952dc1ac5fef89/docs/doc-source/figs/QM_unified_schema_wrapped_280820.gv.png
   :alt: QM Platform schematic

   QM Platform schematic



Citation
--------

If you make use of PEMtk in your research, please cite it.

Cite the software directly via the Github repository for the software - use the "Cite this repository" link in Github, or use the included `CITATION.bib` file, which includes::

  @software{hockett2021PEMtkGithub,
    title = {Photoelectron Metrology Toolkit (PEMtk) Github Repository},
    author = {Hockett, Paul},
    year = {2024},
    url = {https://github.com/phockett/PEMtk},
    urldate = {2022-02-18},
    DOI={10.5281/zenodo.10882996},
    publisher={Github},
    abstract = {Quantum Metrology with Photoelectrons platform data \& analysis layer - the unifiying layer (glue) for the platform. Main capabilities are development of fitting/retrieving continuum wavefunctions from experimental data; handling multi-dimensional datasets; facilitating comparison of ab initio results with experimental data.},
    keywords = {Repo,Software},
    commit = {788329b82911b2a0690323c64116aa6d19537ecc},
  }

(For specific releases and commits, see https://github.com/phockett/PEMtk/releases and https://github.com/phockett/PEMtk/commits/master/, and the archived versions can also be found on Zenodo: |zenodo|.)

... or the book `Quantum Metrology with Photoelectrons Vol. 3 (2023) <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/intro.html>`__ and/or the manuscript `Topical Review: Extracting Molecular Frame Photoionization Dynamics from Experimental Data (2023) <https://www.authorea.com/users/71114/articles/447808-extracting-molecular-frame-photoionization-dynamics-from-experimental-data>`__ paper. Both are available in HTML versions with interactive figures (as linked above), and other forms online, and include discussion and numerical demos using the software (release v0.0.1-dev); see the included `CITATION.bib <https://github.com/phockett/PEMtk/blob/master/CITATION.bib>`__ for additional citation details and options.::

  @book{hockett2023QuantumMetrologyPhotoelectronsIOP,
    title = {Quantum Metrology with Photoelectrons, Volume 3: Analysis Methodologies},
    author = {Hockett, Paul and Makhija, Varun},
    year = {2023},
    month = dec,
    publisher = {IOP Publishing},
    doi = {10.1088/978-0-7503-5022-8},
    url = {https://iopscience.iop.org/book/mono/978-0-7503-5022-8},
    isbn = {978-0-7503-5022-8},
  }

  @article{hockett2023TopicalReviewExtracting,
    title = {Topical Review: Extracting Molecular Frame Photoionization Dynamics from Experimental Data},
    author = {Hockett, Paul and Makhija, Varun},
    year = {2023},
    month = may,
    journal = {Journal of Physics B: Atomic, Molecular and Optical Physics},
    volume = {56},
    number = {11},
    eprint = {2209.04301},
    pages = {112001},
    publisher = {IOP Publishing},
    issn = {0953-4075},
    doi = {10.1088/1361-6455/acd03e},
    url = {https://dx.doi.org/10.1088/1361-6455/acd03e},
  }

(Citation styles for software `from StackExchange <https://academia.stackexchange.com/questions/14010/how-do-you-cite-a-github-repository>`_.)


.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10882996.svg
    :alt: Zenodo archive
    :scale: 100%
    :target: https://doi.org/10.5281/zenodo.10882996
