import pytest
import pyfesom as pf
#from netCDF4 import Dataset
#from collections import OrderedDict
import os

meshpath  ='../FESOM-data/pi-grid/'
mesh = pf.load_mesh(meshpath)

def test_mesh():
    assert mesh.x2.shape[0] == 3140
    assert mesh.y2.shape[0] == 3140
    assert mesh.x2.mean().round(3) == 5.267
    assert mesh.y2.mean().round(3) == 9.729
    assert mesh.elem.shape[0] == 5839
    assert mesh.elem.mean().round(0) == 1581


