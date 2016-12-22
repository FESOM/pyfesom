# pyfesom
# -------
# python tools for FESOM ocean model
#
# Original code by Dmitry Sidorenko, 2013
#
# Copyright (c) 2016 FESOM team

import sys

if sys.version_info < (3,3):
    raise ImportError(
"""
pyfesom does not support Python 2.6, 2.7, 3.0, 3.1, or 3.2.
""")

from .load_mesh_data import *
from .regriding import *
from .woa2005 import *
from .ut import *
from .climatology import *


__author__   = 'FESOM team'
__license__  = 'MIT'
__version__  = '0.1'







