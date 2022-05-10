PEMtk: the Photoelectron Metrology Toolkit
==========================================

(aka Quantum Metrology with Photoelectrons platform data & analysis
layer)

The unifiying layer (glue) for the platform - more details to follow soon.

For now, see

* `Ongoing docs on ReadTheDocs <https://pemtk.readthedocs.io/en/latest/index.html>`__.
* `ePSdata for general aims & motivation <https://phockett.github.io/ePSdata/about.html#Motivation>`__.


.. Local fig: .. figure:: ./docs/doc-source/figs/QM_unified_schema_wrapped_280820.gv.png
   Use GH version via full URL instead for consistency on RTD.

.. figure:: https://raw.githubusercontent.com/phockett/PEMtk/4eec9217203bfd1aee13bd8b64952dc1ac5fef89/docs/doc-source/figs/QM_unified_schema_wrapped_280820.gv.png
   :alt: QM Platform schematic

   QM Platform schematic


Currently implemented
---------------------

- `Symmetrized harmonics generation functions (March 2022)<https://pemtk.readthedocs.io/en/latest/sym/pemtk_symHarm_demo_160322_tidy.html>`_
- [Fitting analysis functions and plotters (Dec. 2021)](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_multiproc_class_analysis_141121-tidy.html).
- Basic data model and [fitting class (20/05/21) + docs](https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html).
- Code framework and notes (August 2020).


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

* The repo can be passed directly to pip, e.g. `pip install git+https://github.com/phockett/PEMtk.git`, see `notes in the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#git>`_.
* Note that `pip -e` is for 'editable', and requires the source dir to remain, but the installation is also editable, `see notes here <https://stackoverflow.com/questions/41535915/python-pip-install-from-local-dir>`_. Alternatively, use `pip install <path_to_local_pkg>`.


Roadmap
-------

- Further fitting methodology & code developments.
- Integrated data handling class(es).
- Image processing.
- Image simulation.
- More versatile plotting routines.
- Interfaces for various experimental platforms.
