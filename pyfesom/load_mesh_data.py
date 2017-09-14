# This file is part of pyfesom
#
################################################################################
#
# Original code by Dmitry Sidorenko, 2013
#
# Modifications:
#   Nikolay Koldunov, 2016
#          - optimisation of reading ASCII fies (add pandas dependency)
#          - move loading and processing of the mesh to the mesh class itself
#
################################################################################

import pandas as pd
import numpy as np
from netCDF4 import Dataset
from . import ut
from scipy.io import netcdf # replace it with netCDF4
import os
import logging
import time
try:
    import cPickle as pickle
except:
    import pickle
import pickle
import sys
import joblib

if (sys.version_info > (3, 0)):
# Python 3 code in this block
    python_version = 3 
else:
# Python 2 code in this block
    python_version = 2


def load_mesh(path, abg = [50, 15, -90], 
              get3d = True, usepickle = True,
              usejoblib = False ):
    ''' Loads FESOM mesh 

    Parameters
    ----------
    path : str
        Path to the directory with mesh files
    abg : list
        alpha, beta and gamma Euler angles. Default [50, 15, -90]
    get3d : bool
        do we load complete 3d mesh or only 2d nodes.

    Returns
    -------
    mesh : object
        fesom_mesh object
    '''
    path=os.path.abspath(path)
    if (usepickle==True) and (usejoblib==True):
        raise ValueError("Both `usepickle` and `usejoblib` set to True, select only one")

    if usepickle:    
        if python_version==3:
            pickle_file = os.path.join(path,'pickle_mesh_py3')
        elif python_version == 2:
            pickle_file = os.path.join(path,'pickle_mesh')
        print(pickle_file)
        print(python_version)

    if usejoblib:
        joblib_file = os.path.join(path, 'joblib_mesh')

    if usepickle and (os.path.isfile(pickle_file)):
        print("The usepickle == True)")
        print("The pickle file for python {} exists.".format(str(python_version)))
        print("The mesh will be loaded from {}".format(pickle_file))

        ifile = open(pickle_file, 'rb')
        mesh = pickle.load(ifile)
        ifile.close()
        return mesh

    elif (usepickle==True) and (os.path.isfile(pickle_file)==False):
        print('The usepickle == True')
        print('The pickle file for python {} DO NOT exists'.format(str(python_version)))
        print('The mesh will be saved to {}'.format(pickle_file))

        mesh = fesom_mesh(path=path, abg=abg, get3d=get3d)
        logging.info('Use pickle to save the mesh information')
        print('Save mesh to binary format')
        outfile = open(pickle_file, 'wb')
        pickle.dump(mesh, outfile)
        outfile.close()
        return mesh

    elif (usepickle==False) and (usejoblib==False):
        mesh = fesom_mesh(path=path, abg=abg, get3d=get3d)
        return mesh

    if (usejoblib==True) and (os.path.isfile(joblib_file)):
        print("The usejoblib == True)")
        print("The joblib file for python {} exists.".format(str(python_version)))
        print("The mesh will be loaded from {}".format(joblib_file))

        mesh = joblib.load(joblib_file)
        return mesh

    elif (usejoblib==True) and (os.path.isfile(joblib_file)==False):
        print('The usejoblib == True')
        print('The joblib file for python {} DO NOT exists'.format(str(python_version)))
        print('The mesh will be saved to {}'.format(joblib_file))

        mesh = fesom_mesh(path=path, abg=abg, get3d=get3d)
        logging.info('Use joblib to save the mesh information')
        print('Save mesh to binary format')
        joblib.dump(mesh, joblib_file)
        
        return mesh


class fesom_mesh(object):
    """ Creates instance of the FESOM mesh.

    This class creates instance that contain information
    about FESOM mesh. At present the class works with 
    ASCII representation of the FESOM grid, but should be extended 
    to be able to read also netCDF version (probably UGRID convention).

    Minimum requirement is to provide the path to the directory, 
    where following files should be located (not nessesarely all of them will
    be used):

    - nod2d.out
    - nod3d.out
    - elem2d.out
    - aux3d.out

    Parameters
    ----------
    path : str
        Path to the directory with mesh files

    abg : list
        alpha, beta and gamma Euler angles. Default [50, 15, -90]

    get3d : bool
        do we load complete 3d mesh or only 2d nodes. 

    Attributes
    ----------
    path : str
        Path to the directory with mesh files
    x2 : array
        x position (lon) of the surface node
    y2 : array
        y position (lat) of the surface node
    n2d : int
        number of 2d nodes
    e2d : int
        number of 2d elements (triangles)
    n3d : int 
        number of 3d nodes
    nlev : int
        number of vertical levels
    zlevs : array
        array of vertical level depths
    voltri : array
        array with 2d volume of triangles
    alpha : float
        Euler angle alpha
    beta : float
        Euler angle beta
    gamma : float
        Euler angle gamma

    Returns
    -------
    mesh : object
        fesom_mesh object
"""
    def __init__(self, path, abg = [50, 15, -90], get3d = True):
        self.path=os.path.abspath(path)

        if not os.path.exists(self.path):
            raise IOError("The path \"{}\" does not exists".format(self.path))


        self.alpha=abg[0]
        self.beta=abg[1]
        self.gamma=abg[2]

        self.nod2dfile=os.path.join(self.path,'nod2d.out')
        self.elm2dfile=os.path.join(self.path,'elem2d.out')
        self.aux3dfile=os.path.join(self.path,'aux3d.out')
        self.nod3dfile=os.path.join(self.path,'nod3d.out')

        self.n3d=0
        self.e2d=0
        self.nlev=0
        self.zlevs= []
        # self.x2= []
        # self.y2= []
        #self.elem= []
        self.n32=[]
        #self.no_cyclic_elem=[]
        self.topo=[]
        self.voltri=[]

        logging.info('load 2d part of the grid')
        start = time.clock()
        self.read2d()
        end = time.clock()
        print('Load 2d part of the grid in {} second(s)'.format(str(int(end-start))))

        if get3d:
            logging.info('load 3d part of the grid')
            start = time.clock()
            self.read3d()
            end = time.clock()
            print('Load 3d part of the grid in {} seconds'.format(str(int(end-start))))
        


    def read2d(self):
        ''' Reads only surface part of the mesh. 
        Useful if your mesh is large and you want to visualize only surface.
        '''
        file_content = pd.read_csv(self.nod2dfile, delim_whitespace=True, skiprows=1, \
                                      names=['node_number','x','y','flag'] )
        self.x2     =file_content.x.values
        self.y2     =file_content.y.values
        self.ind2d  =file_content.flag.values
        self.n2d=len(self.x2)

        file_content = pd.read_csv(self.elm2dfile, delim_whitespace=True, skiprows=1, \
                                              names=['first_elem','second_elem','third_elem'])
        self.elem=file_content.values-1
        self.e2d=np.shape(self.elem)[0]

                ########################################### 
        #here we compute the volumes of the triangles
        #this should be moved into fesom generan mesh output netcdf file
        #
        r_earth=6371000.0
        rad=np.pi/180
        edx=self.x2[self.elem]
        edy=self.y2[self.elem]
        ed=np.array([edx, edy])

        jacobian2D=ed[:, :, 1]-ed[:, :, 0]
        jacobian2D=np.array([jacobian2D, ed[:, :, 2]-ed[:, :, 0]])
        for j in range(2):
            mind = [i for (i, val) in enumerate(jacobian2D[j,0,:]) if val > 355]
            pind = [i for (i, val) in enumerate(jacobian2D[j,0,:]) if val < -355]
            jacobian2D[j,0,mind]=jacobian2D[j,0,mind]-360
            jacobian2D[j,0,pind]=jacobian2D[j,0,pind]+360

        jacobian2D=jacobian2D*r_earth*rad

        for k in range(2):
            jacobian2D[k,0,:]=jacobian2D[k,0,:]*np.cos(edy*rad).mean(axis=1)

        self.voltri = abs(np.linalg.det(np.rollaxis(jacobian2D, 2)))/2.
        
        # compute the 2D lump operator
        cnt=np.array((0,)*self.n2d)
        self.lump2=np.array((0.,)*self.n2d)
        for i in range(3):
            for j in range(self.e2d):
                n=self.elem[j,i]
                #cnt[n]=cnt[n]+1
                self.lump2[n]=self.lump2[n]+self.voltri[j]
        self.lump2=self.lump2/3.

        self.x2, self.y2 = ut.scalar_r2g(self.alpha,self.beta,self.gamma,self.x2, self.y2)
        d=self.x2[self.elem].max(axis=1) - self.x2[self.elem].min(axis=1)
        self.no_cyclic_elem = [i for (i, val) in enumerate(d) if val < 100]

        return self

    def read3d(self):
        '''
        Reads 3d part of the mesh. 
        '''
        self.n3d=int(open(self.nod3dfile).readline().rstrip())
        df = pd.read_csv(self.nod3dfile, delim_whitespace=True, skiprows=1, \
                                    names=['node_number','x','y','z','flag'])
        zcoord = -df.z.values
        self.zlevs = np.unique(zcoord)

        with open(self.aux3dfile) as f:
            self.nlev=int(next(f))
            self.n32=np.array([next(f) for x in \
                            range(self.n2d*self.nlev)]).astype(int).reshape(self.n2d, self.nlev)   
        self.topo=np.zeros(shape=(self.n2d))
        for prof in self.n32:           
            ind_nan = prof[prof>0]
            ind_nan=ind_nan[-1]
            self.topo[prof[0]-1]=zcoord[ind_nan-1]
        
        return self


    def __repr__(self):
        meshinfo = '''
FESOM mesh:
path                  = {}
alpha, beta, gamma    = {}, {}, {}
number of 2d nodes    = {}
number of 2d elements = {}
number of 3d nodes    = {}

        '''.format(self.path,
                   str(self.alpha),
                   str(self.beta),
                   str(self.gamma),
                   str(self.n2d),
                   str(self.e2d),
                   str(self.n3d))
        return meshinfo 
    def __str__(self):
        return __repr__(self)
        

def read_fesom_mesh(path, alpha, beta, gamma, read_diag=True):
    '''
    .. note:: Deprecated in pyfesom 0.1
          `read_fesom_mesh` will be removed in future, it is replaced by
          `load_mesh`. 
          
    ''' 

    mesh=fesom_mesh()
    mesh.path=path
    mesh.alpha=alpha
    mesh.beta=beta
    mesh.gamma=gamma    

    nod2dfile=mesh.path+'nod2d.out'
    elm2dfile=mesh.path+'elem2d.out'
    aux3dfile=mesh.path+'aux3d.out'
    nod3dfile=mesh.path+'nod3d.out'


    file_content = pd.read_csv(nod2dfile, delim_whitespace=True, skiprows=1, \
                                          names=['node_number','x','y','flag'] )
    mesh.x2=file_content.x.values
    mesh.y2=file_content.y.values
    mesh.n2d=len(mesh.x2)
    
    file_content = pd.read_csv(elm2dfile, delim_whitespace=True, skiprows=1, \
                                          names=['first_elem','second_elem','third_elem'])
    mesh.elem=file_content.values-1
    mesh.e2d=np.shape(mesh.elem)[0]

    mesh.n3d=int(open(mesh.path+'nod3d.out').readline().rstrip())
    df = pd.read_csv(nod3dfile, delim_whitespace=True, skiprows=1, \
                                names=['node_number','x','y','z','flag'])
    zcoord = -df.z.values
    mesh.zlevs = np.unique(zcoord)

    with open(aux3dfile) as f:
        mesh.nlev=int(next(f))
        mesh.n32=np.array([next(f) for x in range(mesh.n2d*mesh.nlev)]).astype(int).reshape(mesh.n2d,mesh.nlev)   
    mesh.topo=np.zeros(shape=(mesh.n2d))
    for prof in mesh.n32:           
        ind_nan = prof[prof>0]
        ind_nan=ind_nan[-1]
        mesh.topo[prof[0]-1]=zcoord[ind_nan-1]
    ########################################### 
    #here we compute the volumes of the triangles
    #this should be moved into fesom generan mesh output netcdf file
    #
    r_earth=6371000.0
    rad=np.pi/180
    edx=mesh.x2[mesh.elem]
    edy=mesh.y2[mesh.elem]
    ed=np.array([edx, edy])

    jacobian2D=ed[:, :, 1]-ed[:, :, 0]
    jacobian2D=np.array([jacobian2D, ed[:, :, 2]-ed[:, :, 0]])
    for j in range(2):
        mind = [i for (i, val) in enumerate(jacobian2D[j,0,:]) if val > 355]
        pind = [i for (i, val) in enumerate(jacobian2D[j,0,:]) if val < -355]
        jacobian2D[j,0,mind]=jacobian2D[j,0,mind]-360
        jacobian2D[j,0,pind]=jacobian2D[j,0,pind]+360

    jacobian2D=jacobian2D*r_earth*rad

    for k in range(2):
        jacobian2D[k,0,:]=jacobian2D[k,0,:]*np.cos(edy*rad).mean(axis=1)

    mesh.voltri = abs(np.linalg.det(np.rollaxis(jacobian2D, 2)))/2.
    
    # compute the 2D lump operator
    cnt=np.array((0,)*mesh.n2d)
    mesh.lump2=np.array((0.,)*mesh.n2d)
    for i in range(3):
        for j in range(mesh.e2d):
            n=mesh.elem[j,i]
            #cnt[n]=cnt[n]+1
            mesh.lump2[n]=mesh.lump2[n]+mesh.voltri[j]
    mesh.lump2=mesh.lump2/3.
    ###########################################
    #here we read the 3D cluster volumes
    if (read_diag):
        f = netcdf.netcdf_file(mesh.path+'mesh.diag.nc', 'r')
        mesh.cluster_vol3=f.variables['cluster_vol'].data
        mesh.cluster_vol2=f.variables['cluster_area'].data
        f.close()
    else:
        mesh.cluster_vol3=0
        mesh.cluster_vol2=0
    #we should rotate the mesh to the geographical coordinates
    (mesh.x2,mesh.y2)=ut.scalar_r2g(mesh.alpha,mesh.beta,mesh.gamma,mesh.x2,mesh.y2)
    d=mesh.x2[mesh.elem].max(axis=1)-mesh.x2[mesh.elem].min(axis=1)
    mesh.no_cyclic_elem = [i for (i, val) in enumerate(d) if val < 100]
    return mesh

def read_fesom_3d(str_id, months, years, mesh, result_path, runid, ext, how='mean'): 
    '''
    .. note:: Deprecated in pyfesom 0.1
          `read_fesom_3d` will be removed in future, it is replaced by
          `get_data`. 
          
    '''
    str_id      =str_id
    ext         =ext
    runid       =runid
    years       =years
    months      =months
    result_path =result_path

    y           =years[0]
    data3       =np.zeros(shape=(mesh.n3d))
    while y<=years[1]:
        print(['reading year '+str(y)+':'])
        ncfile =result_path+runid+'.'+str(y)+ext
        f = Dataset(ncfile, 'r')
        if how=='mean':
            data3 = data3+f.variables[str_id][months,:].mean(axis=0)
        elif how=='max':
            data3 = data3+f.variables[str_id][months,:].max(axis=0)
        f.close()
        y=y+1
    data3=data3/(years[1]-years[0]+1)
    return data3

def read_fesom_2d(str_id, months, years, mesh, result_path, runid, ext, how='mean'): 
    '''
    .. note:: Deprecated in pyfesom 0.1
          `read_fesom_2d` will be removed in future, it is replaced by
          `get_data`. 
    '''   
    y=years[0]
    data2=np.zeros(shape=(mesh.n2d))
    while y<=years[1]:
        print(['reading year '+str(y)+':'])
        ncfile =result_path+runid+'.'+str(y)+ext
        f = Dataset(ncfile, 'r')
        if how=='mean':
            data2 = data2+f.variables[str_id][months,0:mesh.n2d].mean(axis=0)
        elif how=='max':
            data2 = data2+f.variables[str_id][months,0:mesh.n2d].max(axis=0)
        f.close()
        y=y+1
    data2=data2/(years[1]-years[0]+1)
    return data2

def ind_for_depth(depth, mesh):
    '''
    Return indexes that belong to certain depth.

    Parameters
    ----------
    depth : float
        desired depth. Note thre will be no interpolation, the model level
        that is closest to desired depth will be selected.
    mesh : object
        FESOM mesh object

    Returns
    -------
    ind_depth : 1d array
        vector with the size equal to the size of the surface nodes with index values where
        we have data values and missing values where we don't have data values.
    ind_noempty : 1d array
        vector with indexes of the `ind_depth` that have data values.
    ind_empty : 1d array
        vector with indexes of the `ind_depth` that do not have data values.
    '''
     
    # Find indexes of the model depth that are closest to the required depth.
    dind=(abs(mesh.zlevs-depth)).argmin()
    # select data from the level and find indexes with values and with nans
    ind_depth=mesh.n32[:,dind]-1
    ind_noempty=np.where(ind_depth>=0)[0]
    ind_empty=np.where(ind_depth<0)[0]
    return ind_depth, ind_noempty, ind_empty

def fesom2depth(depth, data3, mesh, verbose=True):
    '''
    Return 2d array of the 2d mesh shape with data
    from the model level closest to the desired depth. 
    There is no interpolation to the depth.

    Parameters
    ----------
    depth : int
        desired depth
    data3 : array
        complete 3d data (vector) for one timestep
    mesh  : fesom_mesh object
        mesh representation
    verbose : bool
        flag to turn off information about which level will be used.

    Returns
    -------
    data2 : array
        2d array (actually vector) with data from the desired level.


    ''' 

    # create 2d field with the 2d mesh size
    data2=np.zeros(shape=(mesh.n2d))
    dind=(abs(mesh.zlevs-depth)).argmin()
    
    ind_depth, ind_noempty, ind_empty = ind_for_depth(depth, mesh)
    
    # fill in the output array 
    data2[ind_noempty]=data3[ind_depth[ind_noempty]]
    data2[ind_empty]=np.nan
    if verbose:
        print("For depth {} model level {} will be used".format(str(depth),str(mesh.zlevs[dind])))
    return data2



def get_data(data, mesh, depth = 0, verbose=True):
    '''
    Show data from the model level that is closest to the
    desired depth. 

    Parameters
    ----------
    data : array
        complete 3d data for one timestep
    mesh : fesom_mesh object
        mesh representation
    depth : int
        desired depth
    verbose : bool
        flag to turn off information about which level will be used.


    Returns
    -------
    level_data : array
        2d array (actually vector) with data from model level that is closest to the desired depth.
    elem_no_nan : array
        array with triangles (defined as triplets of node indexes) with
        not NaN elements. 

    '''
    elem2=mesh.elem[mesh.no_cyclic_elem,:]
    level_data = fesom2depth(depth, data ,mesh, verbose)
    #The data2[elem2] creates 3d array where every raw is three
    #values of the parameter on the verticies of the triangle.
    d=level_data[elem2].mean(axis=1)
    #k = [i for (i, val) in enumerate(d) if not np.isnan(val)]
    #elem2=elem2[k,:]
    no_nan_triangles = np.invert(np.isnan(d))
    elem_no_nan = elem2[no_nan_triangles,:]

    return level_data, elem_no_nan

def get_layer_mean(data, depth, mesh, timeslice=None):
    '''
    Return mean over the model depth that is closest to specified depth.
    '''

    ind_depth, ind_noempty, ind_empty = ind_for_depth(depth, mesh)
    data_mean=np.zeros(shape=(mesh.n2d))
    if timeslice is None:
        data_mean[ind_noempty]=data[:,ind_depth[ind_noempty]].mean(axis=0)
    else:
        data_mean[ind_noempty]=data[timeslice,ind_depth[ind_noempty]].mean(axis=0)

    data_mean[ind_empty] = np.nan
    #np.ma.masked_equal(data_mean,-9999)

    elem2=mesh.elem[mesh.no_cyclic_elem,:]
    #The data2[elem2] creates 3d array where every raw is three
    #values of the parameter on the verticies of the triangle.
    d=data_mean[elem2].mean(axis=1)
    #k = [i for (i, val) in enumerate(d) if not np.isnan(val)]
    #elem2=elem2[k,:]
    no_nan_triangles = np.invert(np.isnan(d))
    elem_no_nan = elem2[no_nan_triangles,:]
    


    return data_mean, elem_no_nan

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
    ind_depth, ind_noempty, ind_empty = ind_for_depth(depth, mesh)
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



