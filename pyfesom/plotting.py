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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cmocean import cm as cmo
import palettable
import matplotlib as mpl
from .transect import *


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
                  levels = None,
                  cmap = cm.Spectral_r,
                  label = '$^{\circ}$C',
                  title = ''):

    lonlat = transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=npoints)
    nodes  = transect_get_nodes(lonlat, mesh)
    dist   = transect_get_distance(lonlat)
    profile = transect_get_profile(nodes, mesh)
#     data3d =  fl.variables['salt'][0,:]
    transect_data = transect_get_data(data3d, profile)
    
#     plt.figure(figsize=(10,5))
    depth_index=(abs(mesh.zlevs-maxdepth)).argmin()
    plt.subplot(1,1,1)
    trsplt = plt.contourf( dist, mesh.zlevs[:depth_index], transect_data[:,:depth_index].T, levels = levels, 
                 cmap = cmap, 
                 extend='both')
    cb = plt.colorbar()
    cb.set_label(label)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel('km')
    plt.ylabel('m')
    return trsplt