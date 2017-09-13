pyfesom
=======

Python library and collection of tools for basic handling of `FESOM <http://www.fesom.de/>`_  ocean model output.

:ref:`tools` are python scripts with command line interfaces that are used for quick actions with FESOM model output. For example::

    python showme.py /path/to/mesh /path/to/file.nc salt

will produce a map with global spatial distribution of salinity at the surface during the first time step.  

Library is a python library that contains functions for working with FESOM mesh and data. For example loading FESOM mesh can be done as simple as::

    import pyfesom as pf 
    meshpath  ='/path/to/mesh/'
    mesh = pf.load_mesh(meshpath)

Examples of tools are :ref:`showme` for quick visualization of FESOM data and :ref:`scalar2geo` for interpolation to regular lon/lat grid.

.. toctree::
   :maxdepth: 2
   :caption: Content:

   installation
   tools
   library
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
