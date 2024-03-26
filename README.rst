PEMtk: the Photoelectron Metrology Toolkit
==========================================

(aka Quantum Metrology with Photoelectrons platform data & analysis
layer: the unifiying layer (glue) for the platform)

Currently, the platform is under development, based on methods developed over the last 10+ years of research in this area. On the theory side, work is based around ePolyScat, and a set of python packages have already been developed (`ePSproc <https://epsproc.readthedocs.io>`__); on the experimental side, the plan is to update existing Matlab codes for Velocity Map Imaging (VMI) experiments and analysis routines (and rewrite/unify in python). The real foundation, and glue, for the platform will be the Photoelectron metrology toolkit (PEMtk), which will provide the unifying data platform, and analysis routines. In the future, it is hoped that this platform will be extended to other theoretical and experimental methods, but continue to provide a useful, unifying, platform.


As of late 2023, full modelling and fitting of AF and MF observables is implemented, for details see:

* `Ongoing PEMtk docs on ReadTheDocs <https://pemtk.readthedocs.io/en/latest/index.html>`__ - includes getting started guide and basic demos.
* `Quantum Metrology with Photoelectrons Vol. 3 <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/intro.html>`__ - open source book covering the theory and numerical implementation in detail, including Jupyter notebooks and full numerical examples of matrix element retrieval for three demonstration cases (N2, OCS, C2H4), using `v0.0.1-dev release <https://github.com/phockett/PEMtk/releases/tag/v0.0.1-dev-QM3-310723>`__ of the codebase.
* `Topical Review: Extracting Molecular Frame Photoionization Dynamics from Experimental Data (2022) <https://www.authorea.com/users/71114/articles/447808-extracting-molecular-frame-photoionization-dynamics-from-experimental-data>`__ - manuscript which covers photoionization theory and fitting methodologies, the latter as somewhat of a manual for PEMtk, including a full numerical case study.
* See also `ePSdata for general aims & motivation <https://phockett.github.io/ePSdata/about.html#Motivation>`__, and computational results for various cases.

.. Local fig: .. figure:: ./docs/doc-source/figs/QM_unified_schema_wrapped_280820.gv.png
   Use GH version via full URL instead for consistency on RTD.

.. figure:: https://raw.githubusercontent.com/phockett/PEMtk/4eec9217203bfd1aee13bd8b64952dc1ac5fef89/docs/doc-source/figs/QM_unified_schema_wrapped_280820.gv.png
   :alt: QM Platform schematic

   QM Platform schematic


Currently implemented (v0.0.1)
------------------------------

- `Symmetrized harmonics generation functions (March 2022) <https://pemtk.readthedocs.io/en/latest/sym/pemtk_symHarm_demo_160322_tidy.html>`__.
- `Fitting analysis functions and plotters (Dec. 2021) <https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy.html>`__.
- Basic data model (from `ePSproc base class <https://epsproc.readthedocs.io/en/latest/demos/ePSproc_class_demo_161020.html>`__ currently) and `fitting class (20/05/21) + docs <https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html>`__.
- Code framework and notes (August 2020).

A full demonstration can be found in `Quantum Metrology with Photoelectrons Vol. 3 <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/intro.html>`__, which made use of the `v0.0.1-dev release of the codebase (Aug 2023) <https://github.com/phockett/PEMtk/releases/tag/v0.0.1-dev-QM3-310723>`__.


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

* The repo can be passed directly to pip, e.g. ``pip install git+https://github.com/phockett/PEMtk.git``, see `notes in the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#git>`_.
* Note that ``pip -e`` is for 'editable', and requires the source dir to remain, but the installation is also editable, `see notes here <https://stackoverflow.com/questions/41535915/python-pip-install-from-local-dir>`_. Alternatively, use ``pip install <path_to_local_pkg>``.

Docker
------

Docker builds are available as follows:

- `Full image build via Dockerhub for the QM3 book (Aug 2023) <https://hub.docker.com/r/epsproc/quantum-met-vol3>__`. For usage details see the `QM3 Docker builds notes <https://github.com/phockett/Quantum-Metrology-with-Photoelectrons-Vol3#docker-builds>`__. (Also available `via Zenodo DOI:10.5281/zenodo.8286020 <https://doi.org/10.5281/zenodo.8286020>`__.)
- Source Dockerfiles for ePSproc + PEMtk builds from `the Open Photoionization Docker Stacks <https://github.com/phockett/open-photoionization-docker-stacks/tree/main/epsproc-pemtk>`__.


Roadmap
-------

- Further fitting methodology & code developments.
- Integrated data handling class(es).
- Image processing (basic inversions, tomography, FT methods etc.).
- Image simulation.
- More versatile plotting routines.
- Interfaces for various experimental platforms.


Citation
--------

If you make use of PEMtk in your research, please cite it.

Cite the software directly via the Github repository for the software - use the "Cite this repository" link in Github, or use the included `CITATION.bib` file, which includes::

  @software{hockett2021PEMtkGithub,
    title = {Photoelectron Metrology Toolkit (PEMtk) Github Repository},
    author = {Hockett, Paul},
    year = {2021},
    url = {https://github.com/phockett/PEMtk},
    urldate = {2022-02-18},
    DOI={10.5281/zenodo.10044679},
    publisher={Github},
    abstract = {Quantum Metrology with Photoelectrons platform data \& analysis layer - the unifiying layer (glue) for the platform. Main capabilities are development of fitting/retrieving continuum wavefunctions from experimental data; handling multi-dimensional datasets; facilitating comparison of ab initio results with experimental data.},
    keywords = {Repo,Software},
    commit = {788329b82911b2a0690323c64116aa6d19537ecc},
  }

(For specific releases and commits, see https://github.com/phockett/PEMtk/releases and https://github.com/phockett/PEMtk/commits/master/.)

... or the book `Quantum Metrology with Photoelectrons Vol. 3 (2023) <https://phockett.github.io/Quantum-Metrology-with-Photoelectrons-Vol3/intro.html>`__ and/or the manuscript `Topical Review: Extracting Molecular Frame Photoionization Dynamics from Experimental Data (2023) <https://www.authorea.com/users/71114/articles/447808-extracting-molecular-frame-photoionization-dynamics-from-experimental-data>`__ paper. Both include discussion and numerical demos using the software (release v0.0.1-dev), and are available in various flavours online, see the included `CITATION.bib` for citation details and options.

(Citation styles for software `from StackExchange <https://academia.stackexchange.com/questions/14010/how-do-you-cite-a-github-repository>`_.)
