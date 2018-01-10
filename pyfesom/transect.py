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
from .ut import tunnel_fast1d, vec_rotate_r2g
import numpy as np


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

def transect_uv(udata3d, vdata3d, wdata3d,
                mesh, lon_start, lat_start, lon_end, lat_end, npoints=30, abg=[50, 15, -90], myangle = 0):
    '''
    Example:
    rot_u, rot_v, dist, profile = pf.transect_uv(fl.variables['u'][0,:], fl.variables['v'][0,:], fl.variables['w'][0,:],
           mesh, lon_start, lat_start, lon_end, lat_end, myangle=0)

    pf.plot_transect(fl.variables['u'][0,:], mesh, lon_start, lat_start , lon_end, lat_end,
                 transect_data=rot_u.T, dist=dist, profile=profile, 
                 levels= np.round(np.linspace(-0.05, 0.05, 42),4), cmap=cmo.balance, maxdepth=600)
    '''
    lonlat = transect_get_lonlat(lon_start, lat_start, lon_end, lat_end, npoints=40)
    nodes  = transect_get_nodes(lonlat, mesh)
    dist   = transect_get_distance(lonlat)
    profile = transect_get_profile(nodes, mesh)
    
    u = udata3d[profile.flatten()]
    v = vdata3d[profile.flatten()]
    w = wdata3d[profile.flatten()]

    u = u.reshape(profile.shape)
    v = v.reshape(profile.shape)
    w = w.reshape(profile.shape)
    w = np.ma.masked_where(profile==-1000, w)
    
    rot_u = []
    rot_v = []
    for i in range(u.shape[1]):
        uu, vv = vec_rotate_r2g(abg[0], abg[1], abg[2], 
                                   mesh.x2[nodes],
                                   mesh.y2[nodes], 
                                   u[:,i],
                                   v[:,i],
                                   flag=1)
        rot_u.append(uu)
        rot_v.append(vv)
    rot_u = np.array(rot_u)
    rot_v = np.array(rot_v)

    if myangle !=0:
        direct = np.rad2deg(np.arctan2(rot_v,rot_u))
        speed_rot = np.hypot(rot_u, rot_v)

        myangle = myangle
        U = speed_rot * np.cos(np.deg2rad(myangle-direct))
        V = speed_rot * np.sin(np.deg2rad(myangle-direct))

        U = np.ma.masked_where(profile.T==-1000, U)
        V = np.ma.masked_where(profile.T==-1000, V)
        rot_u = U
        rot_v = V
    else:
        return rot_u, rot_v, dist, profile
    return rot_u, rot_v, dist, profile
