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
#
################################################################################
from .load_mesh_data import fesom2depth
from .regriding import fesom2regular, create_indexes_and_distances
import numpy as np
import scipy as sc
from numpy import nanmean
from netCDF4 import Dataset
class woa2005:
    """ Class that contains information from the WOA2005 (TS)


     existing instances are: x,y,z,T,S
    example: w=woa2005(woa05_path)"""
    def __init__(self, woa05_path):
        ncfile= Dataset(woa05_path+'t00an1.nc', 'r')
        self.T= np.copy(ncfile.variables['t00an1'][0,:,:,:])
        x=np.copy(ncfile.variables['lon'][:])
        x[x>180]=x[x>180]-360
        ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
        x=np.sort(x)
        self.x=x
        self.y=ncfile.variables['lat'][:]
        self.z=ncfile.variables['depth'][:]
        ncfile.close()
        self.T[:,:,:]=self.T[:,:,ind]
        ncfile=Dataset(woa05_path+'s00an1.nc', 'r')
        self.S=np.copy(ncfile.variables['s00an1'][0,:,:,:])
        self.S[:,:,:]=self.S[:,:,ind]
        ncfile.close()
        self.T[self.T>90]=np.nan
        self.S[self.S>90]=np.nan
        self.Tyz=nanmean(self.T, 2)
        self.Syz=nanmean(self.S, 2)

def fesom_2_woa2005(data, mesh, woa05, verbose=True):
    '''
    Interpolation of fesom data to WOA2005 grid.
    '''
    xx,yy = np.meshgrid(woa05.x, woa05.y)
    zz=np.copy(woa05.T)
    distances, inds = create_indexes_and_distances(mesh, xx, yy,\
                                                k=10, n_jobs=2)
    for dep_ind in range(len(woa05.z)):
        if verbose:
            print('interpolating level: {}'.format(str(dep_ind)))
        wdep=woa05.z[dep_ind]
        dep_up=[z for z in abs(mesh.zlevs) if z<=wdep][-1]
        dep_lo=[z for z in abs(mesh.zlevs) if z>wdep][0]
        i_up=1-abs(wdep-dep_up)/(dep_lo-dep_up)
        i_lo=1-abs(wdep-dep_lo)/(dep_lo-dep_up)
        data2=i_up*fesom2depth(dep_up, data, mesh, verbose=False)
        data2=data2+i_lo*fesom2depth(dep_lo, data, mesh, verbose=False)
        #zz[dep_ind,:,:] = pf.fesom2regular(data2, mesh, xx,yy)
        zz[dep_ind,:,:] = fesom2regular(data2, mesh, xx, yy, distances=distances,\
                               inds=inds)
    zz[np.isnan(woa05.T)]=np.nan
    return xx, yy, zz

def fesom_2_woa2005_part(data, mesh, woa05, levels=None, verbose=True):

    if levels is None:
        levels = woa05.z
    else:
        levels = np.array(levels)
        check = np.in1d(levels, woa05.z)
        if False in check:
            raise ValueError('Not all of the layers that you specify are WOA2005 layers')

    xx,yy = np.meshgrid(woa05.x, woa05.y)
    zz=np.zeros((levels.shape[0],woa05.T.shape[1], woa05.T.shape[2]))
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
        zz[dep_ind,:,:] = fesom2regular(data2, mesh, xx, yy, distances=distances,\
                               inds=inds)
    depth_indexes = [np.where(woa05.z==i)[0][0] for i in levels]
    zz[np.isnan(woa05.T[depth_indexes,:,:])]=np.nan
    return xx, yy, zz

