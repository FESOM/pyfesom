.. _installation:

Installation
============

Now we support only installation from source and recomend to use `conda <https://conda.io/docs/>`_  to install dependencies. You maight succed with usual :code:`pip install` way, but some dependencies like `cartopy` are pretty hard to install without `conda`. The shortest way to succes consist of the following simple steps:

1. Go to `Miniconda <https://conda.io/miniconda.html>`_ website and download Miniconda installation script for your system. We recomend to use python 2.7 version. 

2. Install `Miniconda`. Don't forget to add the path of `Miniconda` instllation to your `$PATH` and relaunch your terminal.

3. Execute the following lines::

    conda config --add channels conda-forge
    conda install pandas netcdf4 cartopy basemap scipy joblib seawater matplotlib click

4. Go to the folder where you want to have pyfesom and execute (you have to have git installed):

    git clone https://github.com/FESOM/pyfesom.git

Now you hopefully have all dependencies in plate and a version of `pyfesom` on your system. At this point you should be able to use :ref:`tools`.

Tools
-----

Pyfesom tools are simple python scripts that are buld with use of the pyfesom library. To run the tool you should usually execute something like this::

    python /path/to/installation/pyfesom/tools/showme.py /path/to/mesh/ /path/to/file.nc 

That's a lot of letters. To make live easier it is recommended for Linux and Mac OS users to create an alias for every tool. For `bash` users, edit your `.bashrc` (or .bash_profile on Mac)::

    alias showme='python /path/to/installation/pyfesom/tools/showme.py'

For `csh` users edit your `.cshrc`::

    alias showme python /path/to/installation/pyfesom/tools/showme.py

Don't forget to `source` your configuration file afterwards.

If you setup an alias as described above the call to the `showme` tool become::

    showme /path/to/mesh/ /path/to/file.nc

or with some options::

    showme -m merc -d 100 -l -6 6 21 -b -100 20 0 65 /path/to/mesh/ /path/to/file.nc

It also make sence to create system variables for paths to meshes.

Library
-------

Since we do not support the standard installation of pyfesom yet, the easiest way to use the library is to add the path with library location to the system path in the beggining of the script or Jupyter notebook::

    import sys
    sys.path.append("/path/to/installation/pyfesom/")
    import pyfesom as pf



