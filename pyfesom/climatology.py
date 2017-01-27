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
from .regriding import fesom2regular, create_indexes_and_distances
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


def fesom_2_clim(data, mesh, climatology, verbose=True):
    '''
    Interpolation of fesom data to grid of the climatology.

    Parameters
    ----------
    data : array
        1d array of FESOM 3d data for one time step
    mesh : mesh object
        FESOM mesh object
    climatology: climatology object
        FESOM climatology object

    Returns
    -------
    xx : 2d array
        longitudes
    yy : 2d array
        latitudes
    out_data : 2d array
       array with data interpolated to climatology level

    '''
    xx,yy = np.meshgrid(climatology.x, climatology.y)
    out_data=np.copy(climatology.T)
    distances, inds = create_indexes_and_distances(mesh, xx, yy,\
                                                k=10, n_jobs=2)
    for dep_ind in range(len(climatology.z)):
        if verbose:
            print('interpolating level: {}'.format(str(dep_ind)))
        wdep=climatology.z[dep_ind]
        dep_up=[z for z in abs(mesh.zlevs) if z<=wdep][-1]
        dep_lo=[z for z in abs(mesh.zlevs) if z>wdep][0]
        i_up=1-abs(wdep-dep_up)/(dep_lo-dep_up)
        i_lo=1-abs(wdep-dep_lo)/(dep_lo-dep_up)
        data2=i_up*fesom2depth(dep_up, data, mesh, verbose=False)
        data2=data2+i_lo*fesom2depth(dep_lo, data, mesh, verbose=False)
        out_data[dep_ind,:,:] = fesom2regular(data2, mesh, xx, yy, distances=distances,\
                               inds=inds)
    out_data[np.isnan(climatology.T)]=np.nan
    return xx, yy, out_data


def fesom_2_clim_onelevel(data, mesh, climatology, levels=None, verbose=True):
    '''
    Interpolation of fesom data to grid of the climatology for set of levels.

    Parameters
    ----------
    data : array
        1d array of FESOM 3d data for one time step
    mesh : mesh object
        FESOM mesh object
    climatology: climatology object
        FESOM climatology object
    levels : list like
        list of levels to interpolate. At present you can use only
        standard levels of WOA.

    Returns
    -------
    xx : 2d array
        longitudes
    yy : 2d array
        latitudes
    out_data : 2d array
       array with data interpolated to climatology level

    '''    

    if levels is None:
        levels = climatology.z
    else:
        levels = np.array(levels)
        check = np.in1d(levels, climatology.z)
        if False in check:
            raise ValueError('Not all of the layers that you specify are WOA2005 layers')

    xx,yy = np.meshgrid(climatology.x, climatology.y)
    out_data=np.zeros((levels.shape[0],climatology.T.shape[1], climatology.T.shape[2]))
    distances, inds = create_indexes_and_distances(mesh, xx, yy,\
                                                k=10, n_jobs=2)
    for dep_ind in range(len(levels)):
        if verbose:
            print('interpolating level: {}'.format(str(dep_ind)))
        wdep=levels[dep_ind]
        dep_up=[z for z in abs(mesh.zlevs) if z<=wdep][-1]
        dep_lo=[z for z in abs(mesh.zlevs) if z>wdep][0]
        i_up=1-abs(wdep-dep_up)/(dep_lo-dep_up)
        i_lo=1-abs(wdep-dep_lo)/(dep_lo-dep_up)
        data2=i_up*fesom2depth(dep_up, data, mesh, verbose=False)
        data2=data2+i_lo*fesom2depth(dep_lo, data, mesh, verbose=False)
        #zz[dep_ind,:,:] = pf.fesom2regular(data2, mesh, xx,yy)
        out_data[dep_ind,:,:] = fesom2regular(data2, mesh, xx, yy, distances=distances,\
                               inds=inds)
    depth_indexes = [np.where(climatology.z==i)[0][0] for i in levels]
    out_data[np.isnan(climatology.T[depth_indexes,:,:])]=np.nan
    return xx, yy, out_data




