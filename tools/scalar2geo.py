import sys
import os
from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
import pyfesom as pf
import joblib
from joblib import Parallel, delayed
import json
from collections import OrderedDict
import click
from mpl_toolkits.basemap import maskoceans
from scipy.interpolate import griddata
import scipy.spatial.qhull as qhull
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator

@click.command()
@click.argument('meshpath', type=click.Path(exists=True), required=True)
@click.argument('ipath', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('opath', nargs=1, required=False, default='./')
@click.argument('variable', nargs=1, required=False, default='temp')
@click.option('--depths', '-d', default=-1, type=click.STRING,show_default=True,
               help='Depths in meters.')
@click.option('--box', '-b',
              nargs=4,
              type=(click.IntRange(-180, 180),
                    click.IntRange(-180, 180),
                    click.IntRange(-90, 90),
                    click.IntRange(-90, 90)),
             default=(-180,180,-80,90), show_default=True,
             help='Map boundaries in -180 180 -90 90 format.')
@click.option('--res', '-r', nargs=2,
              type=(click.INT, click.INT),
              default=(360, 170), show_default=True,
              help='Number of points along each axis (for lon and  lat).')
@click.option('--influence','-i', default=80000, show_default=True,
              help='Radius of influence for interpolation, in meters.')
@click.option('--interp', type=click.Choice(['nn', 'idist', 'linear', 'cubic']),
              default='nn')
@click.option('--timestep', '-t', default=0, show_default=True,
              help='Timstep from netCDF variable, strats with 0.\
              If -1, all timesteps of the netCDF file will be used.')
@click.option('--abg', nargs=3, type=(click.FLOAT,
                    click.FLOAT,
                    click.FLOAT), default=(50, 15, -90))
def convert(meshpath, ipath, opath, variable, depths, box,
            res, influence, timestep, abg, interp):
    print(ipath)
    mesh = pf.load_mesh(meshpath, abg=abg, usepickle=False, usejoblib=True)
    
    sstep = timestep
    radius_of_influence = influence

    left, right, down, up = box
    lonNumber, latNumber = res

    lonreg = np.linspace(left, right, lonNumber)
    latreg = np.linspace(down, up, latNumber)
    lonreg2, latreg2 = np.meshgrid(lonreg, latreg)
    
    localdir = os.path.join(os.path.dirname(__file__))
    print(localdir)
    with open(localdir+'/CMIP6_Omon.json') as data_file:
        cmore_table = json.load(data_file, object_pairs_hook=OrderedDict)

    with open(localdir+'/CMIP6_SIday.json') as data_file:    
        cmore_table_ice = json.load(data_file, object_pairs_hook=OrderedDict)
    
    depths = np.array(depths.split(' '),dtype='float32')
    if depths[0] == -1:
        dind = range(mesh.zlevs.shape[0])
        realdepth = mesh.zlevs
    else:
        dind = []
        realdepth = []
        for depth in depths:
            ddepth = abs(mesh.zlevs-depth).argmin()
            dind.append(ddepth)
            realdepth.append(mesh.zlevs[ddepth])
    print(dind)
    print(realdepth)
    #realdepth = mesh.zlevs[dind]
    k = 10
    distances, inds = pf.create_indexes_and_distances(mesh, lonreg2, latreg2,\
                                                      k=k, n_jobs=4)
    
    ind_depth_all = []
    ind_noempty_all = []
    ind_empty_all = []

    for i in range(len(mesh.zlevs)):
        ind_depth, ind_noempty, ind_empty = pf.ind_for_depth(mesh.zlevs[i], mesh)
        ind_depth_all.append(ind_depth)
        ind_noempty_all.append(ind_noempty)
        ind_empty_all.append(ind_empty)
    if interp == 'nn':
        topo_interp = pf.fesom2regular(mesh.topo, mesh, lonreg2, latreg2, distances=distances,
                                inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
        k = 1
        distances, inds = pf.create_indexes_and_distances(mesh, lonreg2, latreg2,\
                                                          k=k, n_jobs=4)
        points, qh = None, None
    elif interp == 'idist':
        topo_interp = pf.fesom2regular(mesh.topo, mesh, lonreg2, latreg2, distances=distances,
                                inds=inds, radius_of_influence=radius_of_influence, n_jobs=1, how='idist')
        k = 5
        distances, inds = pf.create_indexes_and_distances(mesh, lonreg2, latreg2,\
                                                          k=k, n_jobs=4)
        points, qh = None, None
    elif interp == 'linear':
        points = np.vstack((mesh.x2, mesh.y2)).T
        qh = qhull.Delaunay(points)
        topo_interp = LinearNDInterpolator(qh, mesh.topo)((lonreg2, latreg2))
        distances, inds = None, None
    elif interp == 'cubic':
        points = np.vstack((mesh.x2, mesh.y2)).T
        qh = qhull.Delaunay(points)
        topo_interp = CloughTocher2DInterpolator(qh, mesh.topo)((lonreg2, latreg2))
        distances, inds = None, None

    mdata = maskoceans(lonreg2, latreg2, topo_interp, resolution = 'h', inlands=False)
    topo = np.ma.masked_where(~mdata.mask, topo_interp)
    ncore=2
    Parallel(n_jobs=ncore, verbose=50)(delayed(scalar2geo)(ifile, opath, variable,
                                       mesh, ind_noempty_all,
                                       ind_empty_all,ind_depth_all, cmore_table, lonreg2, latreg2, 
                                       distances, inds, radius_of_influence, topo, points, interp, qh, timestep, dind, realdepth) for ifile in ipath)
    # scalar2geo(ipath, opath, variable,
    #            mesh, ind_noempty_all,
    #            ind_empty_all,ind_depth_all, cmore_table, lonreg2, latreg2, 
    #            distances, inds, radius_of_influence, topo, points, interp, qh, timestep)

def scalar2geo(ifile, opath, variable,
               mesh, ind_noempty_all,
               ind_empty_all,ind_depth_all, cmore_table, 
               lonreg2, latreg2, distances, inds, radius_of_influence, 
               topo, points, interp, qh, timestep, dind, realdepth):
    print(ifile)
    ext = variable
    #ifile = ipath
    ofile = os.path.join(opath, '{}_{}.nc'.format(os.path.basename(ifile)[:-3], ext))
    
    fl = Dataset(ifile)
    fw = Dataset(ofile, mode='w',data_model='NETCDF4_CLASSIC', )

    fw.createDimension('latitude', lonreg2.shape[0])
    fw.createDimension('longitude', latreg2.shape[1])
    fw.createDimension('time', None)
    fw.createDimension('depth_coord',  len(realdepth) )

    lat = fw.createVariable('latitude', 'd', ('latitude'))
    lat.setncatts(noempty_dict(cmore_table['axis_entry']['latitude']))
    lat[:] = latreg2[:,0].flatten()

    lon = fw.createVariable('longitude', 'd', ('longitude'))
    lon.setncatts(noempty_dict(cmore_table['axis_entry']['longitude']))
    lon[:] = lonreg2[0,:].flatten()

    depth = fw.createVariable('depth_coord','d',('depth_coord'))
    depth.setncatts(noempty_dict(cmore_table['axis_entry']['depth_coord']))
    depth[:] = realdepth[:]

    time = fw.createVariable('time','d',('time'))
    time.setncatts(noempty_dict(cmore_table['axis_entry']['time']))
    time.units = fl.variables['time'].units
    if timestep == -1:
        time[:] = fl.variables['time'][:]
    else:
        time[:] = fl.variables['time'][timestep]

    if fl.variables[variable].shape[1] == mesh.n2d:
        dim3d = False
    elif fl.variables[variable].shape[1] == mesh.n3d:
        dim3d = True
    else:
        raise ValueError('Variable size {} is not equal to number of 2d ({}) or 3d ({}) nodes'.format(fl.variables[variable].shape[1], mesh.n2d, mesh.n3d))
    
    # ii = LinearNDInterpolator(qh, mesh.topo)
    if timestep == -1:
        timesteps = range(fl.variables[variable].shape[0])
    else:
        timesteps = range(timestep, timestep+1)

    if True:
        temp = fw.createVariable(variable,'d',\
                                ('time','depth_coord','latitude','longitude'), \
                                fill_value=-9999, zlib=False, complevel=1)
        
        for ttime in timesteps:

            all_layers = fl.variables[variable][ttime,:]
            level_data=np.zeros(shape=(mesh.n2d))
            inter_data=np.zeros(shape=(len(mesh.zlevs),lonreg2.shape[0], lonreg2.shape[1]))
            # inter_data=np.ma.masked_where(topo[0,:,:].mask, inter_data)
            #for i in range(len(mesh.zlevs)):
            for n, i in enumerate(dind):
                #level_data=np.zeros(shape=(mesh.n2d))
                level_data[ind_noempty_all[i]]=all_layers[ind_depth_all[i][ind_noempty_all[i]]]
                level_data[ind_empty_all[i]] = np.nan
                if interp == 'nn':
                    air_nearest = pf.fesom2regular(level_data, mesh, lonreg2, latreg2, distances=distances,
                                            inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
                elif interp == 'idist':
                    air_nearest = pf.fesom2regular(level_data, mesh, lonreg2, latreg2, distances=distances,
                                            inds=inds, radius_of_influence=radius_of_influence, n_jobs=1, how='idist')
                elif interp == 'linear':
                    air_nearest = LinearNDInterpolator(qh, level_data)((lonreg2, latreg2))
                elif interp == 'cubic':
                    air_nearest =CloughTocher2DInterpolator(qh, level_data)((lonreg2, latreg2))

                air_nearest = np.ma.masked_where(topo.mask, air_nearest)
                temp[ttime,n,:,:] = air_nearest[:,:].filled(-9999)

                print i

    fl.close()
    fw.close()
    
    
    # var2d = 0
    # var3d = 0
    # for varname in out_vars:
    #     if vardir[varname]['dims'] == '2D':
    #         var2d += 1
    #     elif vardir[varname]['dims'] == '3D':
    #         var3d += 1
    # var3d = var3d*len(levels)*fl.variables['time'].shape[0]
    # var2d = var2d*fl.variables['time'].shape[0]
    # progress_total  = var3d+var2d 
    # progress_passed = 0

def noempty_dict(d):
    '''
    Removes keys with empty string values from dictionary.

    Parameters
    ----------
    d : OrderedDict
        input dictionary

    Returns
    -------
    d_out : OrderedDict
        output dict with empty strings removed
    '''
    d_out = OrderedDict()
    for key, value in d.iteritems():
        if value != u'':
            d_out[key]=value
    return d_out

if __name__ == '__main__':
    convert()
