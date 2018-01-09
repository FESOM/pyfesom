# This file is part of pyfesom
#
################################################################################
#
# Original code by Nikolay Koldunov, 2018
#
#
################################################################################


import numpy as np
import pyproj
from .ut import tunnel_fast1d

g = pyproj.Geod(ellps='WGS84')

def transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=30):
    lonlat = g.npts(lon_start, lat_start, lon_end, lat_end, npoints)
    lonlat = np.array(lonlat)
    return lonlat

def transect_get_nodes(lonlat, mesh):
    nodes = []
    for i in range(lonlat.shape[0]):
        nn = tunnel_fast1d(mesh.y2, mesh.x2, lonlat[i,1], lonlat[i,0])
        nodes.append(nn[0])
    return nodes

def transect_get_distance(lonlat):
    (az12, az21, dist) = g.inv(lonlat[:,0][0:-1], lonlat[:,1][0:-1], lonlat[:,0][1:], lonlat[:,1][1:])
    dist = dist.cumsum()/1000
    dist = np.insert(dist, 0, 0)
    return dist

def transect_get_profile(nodes, mesh):
    profile = (mesh.n32-1)[nodes,:]
    return profile

def transect_get_data(data3d, profile):
    transect_data = data3d[profile.flatten()]
    transect_data = transect_data.reshape(profile.shape)
    transect_data = np.ma.masked_where(profile==-1000, transect_data)
    return transect_data


