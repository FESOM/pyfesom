# This file is part of pyfesom
#
################################################################################
#
# Original matlab/python code by Sergey Danilov, Dmitry Sidorenko and Qiang Wang.
# 
# Contributers: Lukrecia Stulic, Nikolay Koldunov
#
# Modifications:
#
################################################################################

import numpy as np
import math as mt
import load_mesh_data


def scalar_r2g(al, be, ga, rlon, rlat):
    '''
    Converts rotated coordinates to geographical coordinates.

    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float 
        beta Euler angle
    ga : float 
        gamma Euler angle
    rlon : array
        1d array of longitudes in rotated coordinates
    rlat : array
        1d araay of latitudes in rotated coordinates

    Returns
    -------
    lon : array 
        1d array of longitudes in geographical coordinates
    lat : array
        1d array of latitudes in geographical coordinates

    '''
   
    rad=mt.pi/180
    al=al*rad
    be=be*rad
    ga=ga*rad
    rotate_matrix=np.zeros(shape=(3,3))
    rotate_matrix[0,0]=np.cos(ga)*np.cos(al)-np.sin(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[0,1]=np.cos(ga)*np.sin(al)+np.sin(ga)*np.cos(be)*np.cos(al)
    rotate_matrix[0,2]=np.sin(ga)*np.sin(be)
    rotate_matrix[1,0]=-np.sin(ga)*np.cos(al)-np.cos(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[1,1]=-np.sin(ga)*np.sin(al)+np.cos(ga)*np.cos(be)*np.cos(al)
    rotate_matrix[1,2]=np.cos(ga)*np.sin(be)
    rotate_matrix[2,0]=np.sin(be)*np.sin(al)
    rotate_matrix[2,1]=-np.sin(be)*np.cos(al)
    rotate_matrix[2,2]=np.cos(be)

    rotate_matrix=np.linalg.pinv(rotate_matrix)
    
    rlat=rlat*rad
    rlon=rlon*rad   

    #Rotated Cartesian coordinates:
    xr=np.cos(rlat)*np.cos(rlon)
    yr=np.cos(rlat)*np.sin(rlon)
    zr=np.sin(rlat) 

    #Geographical Cartesian coordinates:
    xg=rotate_matrix[0,0]*xr + rotate_matrix[0,1]*yr + rotate_matrix[0,2]*zr
    yg=rotate_matrix[1,0]*xr + rotate_matrix[1,1]*yr + rotate_matrix[1,2]*zr
    zg=rotate_matrix[2,0]*xr + rotate_matrix[2,1]*yr + rotate_matrix[2,2]*zr        

    #Geographical coordinates:
    lat = np.arcsin(zg)
    lon=  np.arctan2(yg, xg)
    
    a = np.where((np.abs(xg)+np.abs(yg))==0)
    if a: lon[a]=0
    
    lat = lat/rad
    lon = lon/rad

    return (lon,lat)

def scalar_g2r(al, be, ga, lon, lat):
    '''
    Converts geographical coordinates to rotated coordinates.

    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float 
        beta Euler angle
    ga : float 
        gamma Euler angle
    lon : array 
        1d array of longitudes in geographical coordinates
    lat : array
        1d array of latitudes in geographical coordinates
    
    Returns
    -------
    rlon : array
        1d array of longitudes in rotated coordinates
    rlat : array
        1d araay of latitudes in rotated coordinates
    '''

    rad=mt.pi/180
    al=al*rad
    be=be*rad
    ga=ga*rad

    rotate_matrix=np.zeros(shape=(3,3))
    
    rotate_matrix[0,0]=np.cos(ga)*np.cos(al)-np.sin(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[0,1]=np.cos(ga)*np.sin(al)+np.sin(ga)*np.cos(be)*np.cos(al);
    rotate_matrix[0,2]=np.sin(ga)*np.sin(be)
    rotate_matrix[1,0]=-np.sin(ga)*np.cos(al)-np.cos(ga)*np.cos(be)*np.sin(al);
    rotate_matrix[1,1]=-np.sin(ga)*np.sin(al)+np.cos(ga)*np.cos(be)*np.cos(al);
    rotate_matrix[1,2]=np.cos(ga)*np.sin(be);
    rotate_matrix[2,0]=np.sin(be)*np.sin(al);
    rotate_matrix[2,1]=-np.sin(be)*np.cos(al);
    rotate_matrix[2,2]=np.cos(be);
    
    lat=lat*rad;
    lon=lon*rad;

    # geographical Cartesian coordinates:
    xr=np.cos(lat)*np.cos(lon);
    yr=np.cos(lat)*np.sin(lon);
    zr=np.sin(lat);

    # rotated Cartesian coordinates:
    xg=rotate_matrix[0,0]*xr + rotate_matrix[0,1]*yr + rotate_matrix[0,2]*zr;
    yg=rotate_matrix[1,0]*xr + rotate_matrix[1,1]*yr + rotate_matrix[1,2]*zr;
    zg=rotate_matrix[2,0]*xr + rotate_matrix[2,1]*yr + rotate_matrix[2,2]*zr;

    # rotated coordinates:
    rlat=np.arcsin(zg)
    rlon=np.arctan2(yg, xg)

    a = np.where((np.abs(xg)+np.abs(yg))==0)
    if a: lon[a]=0

    rlat = rlat/rad 
    rlon = rlon/rad
    
    return (rlon, rlat)


def vec_rotate_r2g(al, be, ga, lon, lat, urot, vrot, flag):
    '''
    Rotate vectors from rotated coordinates to geographical coordinates.

    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float 
        beta Euler angle
    ga : float 
        gamma Euler angle
    lon : array
        1d array of longitudes in rotated or geographical coordinates (see flag parameter)
    lat : array
        1d array of latitudes in rotated or geographical coordinates (see flag parameter)
    urot : array
        1d array of u component of the vector in rotated coordinates
    vrot : array
        1d array of v component of the vector in rotated coordinates
    flag : 1 or 0
        flag=1  - lon,lat are in geographical coordinate
        flag=0  - lon,lat are in rotated coordinate
    
    Returns
    -------
    u : array
        1d array of u component of the vector in geographical coordinates
    v : array
        1d array of v component of the vector in geographical coordinates

    '''

#   first get another coordinate
    if (flag==1): 
        (rlon,rlat)=scalar_g2r(al, be, ga, lon, lat)
    else:
        rlon=lon
        rlat=lat
        (lon,lat)=scalar_r2g(al, be, ga, rlon, rlat)
 
#   then proceed...
    rad=mt.pi/180
    al=al*rad
    be=be*rad
    ga=ga*rad

    rotate_matrix=np.zeros(shape=(3,3))
    rotate_matrix[0,0]=np.cos(ga)*np.cos(al)-np.sin(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[0,1]=np.cos(ga)*np.sin(al)+np.sin(ga)*np.cos(be)*np.cos(al)
    rotate_matrix[0,2]=np.sin(ga)*np.sin(be)
    rotate_matrix[1,0]=-np.sin(ga)*np.cos(al)-np.cos(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[1,1]=-np.sin(ga)*np.sin(al)+np.cos(ga)*np.cos(be)*np.cos(al)
    rotate_matrix[1,2]=np.cos(ga)*np.sin(be)
    rotate_matrix[2,0]=np.sin(be)*np.sin(al)
    rotate_matrix[2,1]=-np.sin(be)*np.cos(al)
    rotate_matrix[2,2]=np.cos(be)

    rotate_matrix=np.linalg.pinv(rotate_matrix)
    rlat=rlat*rad
    rlon=rlon*rad	
    lat=lat*rad
    lon=lon*rad

#   vector in rotated Cartesian
    txg=-vrot*np.sin(rlat)*np.cos(rlon)-urot*np.sin(rlon)
    tyg=-vrot*np.sin(rlat)*np.sin(rlon)+urot*np.cos(rlon)
    tzg=vrot*np.cos(rlat)
    
#   vector in geo Cartesian
    txr=rotate_matrix[0,0]*txg + rotate_matrix[0,1]*tyg + rotate_matrix[0,2]*tzg 
    tyr=rotate_matrix[1,0]*txg + rotate_matrix[1,1]*tyg + rotate_matrix[1,2]*tzg 
    tzr=rotate_matrix[2,0]*txg + rotate_matrix[2,1]*tyg + rotate_matrix[2,2]*tzg 
    
#   vector in geo coordinate
    v=-np.sin(lat)*np.cos(lon)*txr - np.sin(lat)*np.sin(lon)*tyr + np.cos(lat)*tzr
    u=-np.sin(lon)*txr + np.cos(lon)*tyr
    
    u=np.array(u)
    v=np.array(v)

    return (u,v)

def cut_region(mesh, box=[13, 30, 53, 66], depth=0):
    '''
    Cut region from the mesh.

    Parameters
    ----------
    mesh : object
        FESOM mesh object
    box : list
        Coordinates of the box in [-180 180 -90 90] format.
        Default set to [13, 30, 53, 66], Baltic Sea.
    depth : float
        depth

    Returns
    -------
    elem_no_nan : array
        elements that belong to the region defined by `box`.
    no_nan_triangles : array
        boolian array of size elem2d with True for elements 
        that belong to the region defines by `box`.
    '''
    depth = 0
    left, right, down, up = box
    ind_depth, ind_noempty, ind_empty = load_mesh_data.ind_for_depth(depth, mesh)
    elem2 = mesh.elem
    xx = mesh.x2[elem2]
    yy = mesh.y2[elem2]
    dd = ind_depth[elem2]

    vvv = ((dd>=0) & (xx>=left) & (xx<=right) & (yy>=down) & (yy<=up))
    indd = np.where((ind_depth>=0) &
                    (mesh.x2>=left) &
                    (mesh.x2<=right) &
                    (mesh.y2>=down) &
                    (mesh.y2<=up))
    ccc = vvv.mean(axis=1)
    ccc[ccc!=1] = np.nan

    no_nan_triangles = np.invert(np.isnan(ccc))
    elem_no_nan = elem2[no_nan_triangles,:]

    return elem_no_nan, no_nan_triangles


