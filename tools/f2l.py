from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
import pyfesom as pf
import joblib
from joblib import Parallel, delayed
import click


@click.command()
@click.argument('meshpath', type=click.Path(exists=True), required=True)
@click.argument('ipath', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('opath', nargs=1, required=False, default='./')
@click.argument('variable', nargs=1, required=False, default='temp')
@click.option('--ncore', '-n', default=2, help = 'Number of cores (parallel processes)', show_default=True)
@click.option('--skip', '-s', is_flag=True,
              help='Skip the calculation if the output file already exist.')
def convert(ipath, opath, variable, ncore, meshpath, skip):
    
    mesh = pf.load_mesh(meshpath, usepickle=False, usejoblib=True)

    ind_depth_all = []
    ind_noempty_all = []
    ind_empty_all = []

    for i in range(len(mesh.zlevs)):
        ind_depth, ind_noempty, ind_empty = pf.ind_for_depth(mesh.zlevs[i], mesh)
        ind_depth_all.append(ind_depth)
        ind_noempty_all.append(ind_noempty)
        ind_empty_all.append(ind_empty)

    Parallel(n_jobs=ncore, verbose=50)(delayed(f2l)(l, opath, variable,
                                                mesh, skip, ind_noempty_all, 
                                                ind_empty_all,ind_depth_all) for l in ipath)



def f2l(ifile, opath, variable, mesh, skip, ind_noempty_all, ind_empty_all,ind_depth_all):
    print(ifile)
    ext = 'levels'
    ofile = os.path.join(opath, '{}_{}.nc'.format(os.path.basename(ifile)[:-3], ext))
    if skip:
        if os.path.isfile(ofile):
            print('File {} exist, --skip flag is present, skipping'.format(ofile))
            return
    else:
        try:
            os.remove(ofile)
        except OSError:
            pass  

    fl = Dataset(ifile)
    fw = Dataset(ofile, mode='w', data_model='NETCDF4' )
    fw.createDimension('time', None)
    fw.createDimension('depth',  len(mesh.zlevs) )
    fw.createDimension('ncells',  mesh.n2d )

    time = fw.createVariable('time', 'd', ('time'))
    var = fw.createVariable(variable, 'f', ('time','depth','ncells'), zlib=False, complevel=1, fill_value=1.e+20)
    depth = fw.createVariable('depth', 'd' ,('depth'))
    # data = fl.variables['temp'][0,:]
    level_data=np.zeros(shape=(len(mesh.zlevs), mesh.n2d))

    llim = fl.variables[variable].shape[0]
    #llim = 12

    time.units = fl.variables['time'].units
    time.calendar = "proleptic_gregorian"
    time.long_name = 'time'
    # substruct 1 second to fix the date
    time[:] = fl.variables['time'][:llim] - 1
    
    depth.units = "m"
    depth.long_name = "depth"
    depth.positive = "down"
    depth.axis = "Z"

    depth[:] = mesh.zlevs
    
    var.description = fl.variables[variable].description
    var.units = fl.variables[variable].units
    var.grid_type = "unstructured" 

    for t in range(llim):
        data = fl.variables[variable][t,:]
        for i in range(len(mesh.zlevs)):
            level_data[i, ind_noempty_all[i]]=data[ind_depth_all[i][ind_noempty_all[i]]]
            level_data[i,ind_empty_all[i]] = 1.e+20
        var[t,:,:] = level_data
        if (t>0) and (t%50 == 0):
            print('sync')
            fl.sync()
    fw.close()
    fl.close()


if __name__ == '__main__':
    convert()


















