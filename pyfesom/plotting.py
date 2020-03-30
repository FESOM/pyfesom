# This file is part of pyfesom
#
################################################################################
#
# Original code by Nikolay Koldunov, 2018
#
#
################################################################################

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except:
    print('cartopy is not avalible, plotting will not work')
from cmocean import cm as cmo
import matplotlib as mpl
from .transect import *
import math


def plot_transect_map(lon_start, lat_start, lon_end, lat_end, 
                      mesh, npoints=30, view = 'w', stock_img=False):
    # plt.figure(figsize=(10,10))
    lonlat = transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=npoints)
    nodes  = transect_get_nodes(lonlat, mesh)
    dist   = transect_get_distance(lonlat)

    if view=='w':
        ax = plt.subplot(111, projection=ccrs.Mercator(central_longitude=0))
        ax.set_extent([180, -180, -80, 90], crs=ccrs.PlateCarree())
    elif view=='np':
        ax = plt.subplot(111, projection=ccrs.NorthPolarStereo(central_longitude=0))
        ax.set_extent([180, -180, 60, 90], crs=ccrs.PlateCarree())
    elif view=='sp':
        ax = plt.subplot(111, projection=ccrs.SouthPolarStereo(central_longitude=0))
        ax.set_extent([180, -180, -90, -50], crs=ccrs.PlateCarree())
    else:
        raise ValueError('The "{}" is not recognized as valid view option.'.format(view))
    
    
    ax.scatter(lonlat[:,0], lonlat[:,1], s=30, c='b',  transform=ccrs.PlateCarree())
    ax.scatter(mesh.x2[nodes], mesh.y2[nodes], s=30, c='r',  transform=ccrs.PlateCarree())
    if stock_img==True:
        ax.stock_img()
    ax.coastlines(resolution='50m')
    return ax

def plot_transect(data3d, mesh, lon_start,
                  lat_start,
                  lon_end,
                  lat_end,
                  npoints=30,
                  maxdepth = 1000,
                  label = '$^{\circ}$C',
                  title = '',
                  levels=None,
                  cmap = cm.Spectral_r,
                  ax = None, 
                  dist = None,
                  profile = None, 
                  ncols = 2,
                  figsize = None, 
                  transect_data = []):

        
    depth_index=(abs(mesh.zlevs-maxdepth)).argmin()
    if not isinstance(data3d, list):
        if ax is None:
            ax = plt.gca()
            oneplot = True
        else:
            oneplot = False
        if ((type(dist) is np.ndarray) and (type(profile) is np.ndarray)):
            if not (type(transect_data) is np.ma.core.MaskedArray):
                transect_data = transect_get_data(data3d, profile)
        elif ((type(dist) is np.ndarray) and (type(transect_data) is np.ndarray)):
            lonlat = transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=npoints)
            
        else:
            lonlat = transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=npoints)
            nodes  = transect_get_nodes(lonlat, mesh)
            dist   = transect_get_distance(lonlat)
            profile = transect_get_profile(nodes, mesh)
                
            if not (type(transect_data) is np.ma.core.MaskedArray):
                transect_data = transect_get_data(data3d, profile)
        image = ax.contourf( dist, mesh.zlevs[:depth_index], transect_data[:,:depth_index].T,
                            levels = levels, cmap = cmap, extend='both')
        ax.invert_yaxis()
        ax.set_title(title)
        ax.set_xlabel('km')
        ax.set_ylabel('m')

        if oneplot:
            cb = plt.colorbar(image)
            cb.set_label(label)


        return image
    else:
            ncols = float(ncols)
            nplots = len(data3d)
            nrows = math.ceil(nplots/ncols)
            ncols = int(ncols)
            nrows = int(nrows)
            nplot = 1

            if not figsize:
                figsize = (8*ncols,2*nrows*ncols)
            fig, ax = plt.subplots(nrows,ncols, figsize=figsize)
            ax = ax.flatten()
            for ind, data in enumerate(data3d):
                if ((type(dist) is np.ndarray) and (type(profile) is np.ndarray)):
                    transect_data = transect_get_data(data, profile)
                else:
                    lonlat = transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=npoints)
                    nodes  = transect_get_nodes(lonlat, mesh)
                    dist   = transect_get_distance(lonlat)
                    profile = transect_get_profile(nodes, mesh)
                    transect_data = transect_get_data(data, profile)

                image = ax[ind].contourf( dist, mesh.zlevs[:depth_index], transect_data[:,:depth_index].T,
                                    levels = levels, cmap = cmap, extend='both')
                ax[ind].invert_yaxis()
                if not isinstance(title, list):
                    ax[ind].set_title(title)
                else:
                    ax[ind].set_title(title[ind])
                ax[ind].set_xlabel('km')
                ax[ind].set_ylabel('m')
                
                cb = fig.colorbar(image, orientation='horizontal', ax=ax[ind], pad=0.11)
                cb.set_label(label)
            for delind in range(ind+1, len(ax)):
                
                fig.delaxes(ax[delind])

            fig.tight_layout()
            
                
