pyfesom
=======

Python library and collection of tools for basic handling of `FESOM <http://www.fesom.de/>`_  ocean model output.

:ref:`tools` are python scripts with command line interfaces that are used for quick actions with FESOM model output. For example::

    python showme.py /path/to/mesh /path/to/file.nc salt

will produce map with global spatial distribution of salinity at the surface during the first time step.  

Library is a classical python library that contains functions for working with FESOM mesh and data. For example loading FESOM mesh can be done as simple as::

    import pyfesom as pf 
    meshpath  ='/path/to/mesh/'
    mesh = pf.load_mesh(meshpath)

Examples of tools are :ref:`showme` for quick visualization of FESOM data and scalar2geo for interpolation

Examples of library usage
=========================

- `Show FESOM mesh`_
- `Plot variable on original grid`_
- `Plot simple diagnostics`_ 
- `Interpolate to regular grid`_
- `Compare to climatology`_


.. _Show FESOM mesh: https://github.com/koldunovn/pyfesom/blob/master/notebooks/show_mesh.ipynb
.. _Plot variable on original grid: https://github.com/koldunovn/pyfesom/blob/master/notebooks/show_variable_on_original_grid.ipynb
.. _Plot simple diagnostics: https://github.com/koldunovn/pyfesom/blob/master/notebooks/plot_simple_diagnostics.ipynb
.. _Interpolate to regular grid: https://github.com/koldunovn/pyfesom/blob/master/notebooks/interpolate_to_regular_grid.ipynb
.. _Compare to climatology: https://github.com/koldunovn/pyfesom/blob/master/notebooks/compare_to_climatology.ipynb


Requirements
============

- numpy
- scipy
- pandas
- netCDF4

.. toctree::
   :maxdepth: 2
   :caption: Content:

   installation
   tools
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
