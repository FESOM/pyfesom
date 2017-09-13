.. _f2l:

f2l
======

Fesom results are saved as 1d arrays, where every value is associated with the mesh node. The connectivity (information on how to connect nodes to triangles) is usually only available for surface nodes. In order to get 2d field at other levels (to, for example, plot it) one usually create a dummy 1d array of the size equal to the number of the surface nodes and fill it with values from the model level where values exist and with NaNs where they don't (e.g. there is no level with certain depth for this surface node). 

The `f2l` do exactly this - it takes 1d field that represents model 3d field and splits it in to levels, so instead of, say [time, all_3d_nodes] the dimentions become [time, level, 2d_surface_nodes]. Naturally number of nodes at the surface will be the largest among all levels, so the resulting file will be much bigger. 

Basic usage
-----------
As minimum you should provide path to the mesh, path to the file, path were the ouptut will be stored and variable name::

    python f2l.py /path/to/mesh/ /path/to/file.nc /path/to/output/ temp

The tool can process several files at once. You just have to specify path to your files with the wildcard like this::

   python f2l.py -n 8 /path/to/mesh/ /path/to/file_year*.nc /path/to/output/ temp

By default only once processor is used, so if you would like files to be processed in parallel, you have to specify the number of parallel processes `-n`, that is usually equals to the number of processors you would like to use.

Usage and options
-----------------

Below you can find complete list of options. You can allways display this list in the terminal by executing::

    python f2l.py --help

::

    Usage: f2l.py [OPTIONS] MESHPATH IPATH... [OPATH] [VARIABLE]

    Options:
    -n, --ncore INTEGER  Number of cores (parallel processes)  [default: 2]
    -s, --skip           Skip the calculation if the output file already exist.
    --help               Show this message and exit.