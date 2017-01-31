# This file is part of pyfesom
#
################################################################################
#
# Original code by Dmitry Sidorenko, 2013
#
# Modifications:
#   Nikolay Koldunov, 2016
#          - change to netCDF4 
#          - change scipy griddata interpolation to KDTree for speed
# TODO
# Add seasonal climatology
################################################################################
from .load_mesh_data import fesom2depth
import numpy as np
import scipy as sc
from numpy import nanmean
from netCDF4 import Dataset
import os

class climatology(object):
    '''
    Class that contains information from ocean 1 degree annual climatology.
    Presently options are WOA2005 and PHC3.0

    Parameters
    ----------
    path : str
        Path to the directory with climatology files
    climname : str
        Name of the climatology ('woa05' or 'phc')

    Returns
    -------
    object with climatology fields

    Attributes
    ----------
    x : 2d array
        longitudes
    y : 2d array
        latitudes
    T : 3d array
        temperatures
    S : 3d array
        salinity
    z : 1d array
        depths
    Tyz : 2d array
        zonal mean of temperature
    Syz : 3d array
        zonal mean of salinity

    '''
    def __init__(self, path, climname='woa05'):
        if climname=='woa05':
            ncfile = Dataset(os.path.join(path, 'woa2005TS.nc'))
            self.T = np.copy(ncfile.variables['t00an1'][0,:,:,:])
            x=np.copy(ncfile.variables['lon'][:])
            x[x>180]=x[x>180]-360
            ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
            x=np.sort(x)
            self.x=x
            self.y=ncfile.variables['lat'][:]
            self.z=ncfile.variables['depth'][:]
            self.T[:,:,:]=self.T[:,:,ind]
            self.S=np.copy(ncfile.variables['s00an1'][0,:,:,:])
            self.S[:,:,:]=self.S[:,:,ind]
            ncfile.close()
            self.Tyz=nanmean(self.T, 2)
            self.Syz=nanmean(self.S, 2)
            self.T = np.ma.masked_greater(self.T,1000)
            self.S = np.ma.masked_greater(self.S,1000)

        if climname=='phc':
            ncfile = Dataset(os.path.join(path, 'phc3.0_annual.nc'))
            self.T = np.copy(ncfile.variables['temp'][:,:,:])
            x=np.copy(ncfile.variables['lon'][:])
            x[x>180]=x[x>180]-360
            ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
            x=np.sort(x)
            self.x=x
            self.y=ncfile.variables['lat'][:]
            self.z=ncfile.variables['depth'][:]
            self.T[:,:,:]=self.T[:,:,ind]
            self.S=np.copy(ncfile.variables['salt'][:,:,:])
            self.S[:,:,:]=self.S[:,:,ind]
            ncfile.close()
            self.Tyz=nanmean(self.T, 2)
            self.Syz=nanmean(self.S, 2)

        if climname=='gdem':
            ncfile = Dataset(os.path.join(path, 'gdemv3s_tm.nc'))
            self.T = np.copy(ncfile.variables['water_temp'][0,:,:,:])
            x=np.copy(ncfile.variables['lon'][:])
            x[x>180]=x[x>180]-360
            ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
            x=np.sort(x)
            self.x=x
            self.y=ncfile.variables['lat'][:]
            self.z=ncfile.variables['depth'][:]
            self.T[:,:,:]=self.T[:,:,ind]
            self.S=np.copy(ncfile.variables['salinity'][0,:,:,:])
            self.S[:,:,:]=self.S[:,:,ind]
            ncfile.close()
            self.Tyz=nanmean(self.T, 2)
            self.Syz=nanmean(self.S, 2)
            self.T = np.ma.masked_less(self.T,-1000)
            self.S = np.ma.masked_less(self.S,-1000)


