# This file is part of fesom_viz
#
################################################################################
# 
# Interpolates FESOM data to regular grid. The output is netCDF4 files
# with CMOR complient attributes.   
# 
# Original code by Nikolay Koldunov, 2016
# 
# TODO:
#     Add possibility to define curvilinear grid.
#     Add ESMPy as regridding tool.
#     Find official CMOR descriptions for some of the variables 
# Modifications:
#
################################################################################

import sys
import ConfigParser
from netCDF4 import Dataset, num2date

import numpy as np
import json
from collections import OrderedDict
import os
import datetime
# Read configuration file
try:
    config_file = sys.argv[1]
except:
    print "You have to provide configuration file. Example config located in \"./configs/fesom2geo_example\""
    sys.exit(1)
import multiprocessing as mp

config = ConfigParser.RawConfigParser()
config.read(config_file)

# There is an option to provide path to the pyfesom folder
pfpath = config.get('main', 'pfpath')

sys.path.append(pfpath)
import pyfesom as pf

# Read options from configuration file. See fesom2geo_example for
# explination and possible values
left_lon       = config.getfloat('main', 'left_lon')
right_lon      = config.getfloat('main', 'right_lon')
number_of_lons = config.getint('main', 'number_of_lons')

lower_lat      = config.getfloat('main', 'lower_lat')
upper_lat      = config.getfloat('main', 'upper_lat')
number_of_lats = config.getint('main', 'number_of_lats')


meshpath       = config.get('main', 'meshpath')
path_to_data   = config.get('main', 'path_to_data')
path_to_output = config.get('main', 'path_to_output')

zlib           = config.getboolean('main', 'zlib')
radius_of_influence = config.getint('main', 'radius_of_influence')
k              = config.getint('main', 'neighboring_points')

out_vars       = config.get('main', 'out_vars').split(',')
out_vars = [w.strip() for w in out_vars]
print('='*50)
print("Variables that will be converted: {}".format(out_vars))

levels = np.asarray(config.get('main', 'levels').split(','), dtype='float')


angles_for_mesh     = list(map(int,config.get('main','angles_for_mesh').split(',')))
angles_for_rotation = list(map(int,config.get('main','angles_for_rotation').split(',')))

start_year     = config.getint('main','start_year')
end_year       = config.getint('main','end_year')

ifile_template = config.get('main','ifile_template')
ifile_template_ice = config.get('main','ifile_template_ice')
ofile_template =config.get('main','ofile_template')

distribute_timesteps = config.getboolean('main', 'distribute_timesteps')

# Generate regular grid
lon = np.linspace(left_lon, right_lon, number_of_lons)
lat = np.linspace(lower_lat, upper_lat, number_of_lats)
lons, lats = np.meshgrid(lon,lat)

# read the FESOM mesh
print('='*50)
mesh = pf.load_mesh(meshpath,abg=angles_for_mesh, get3d=True, usepickle=True)

# Open CMOR variable descriptions
with open('CMIP6_Omon.json') as data_file:    
    cmore_table = json.load(data_file, object_pairs_hook=OrderedDict)

with open('CMIP6_SIday.json') as data_file:    
    cmore_table_ice = json.load(data_file, object_pairs_hook=OrderedDict)

# Add some variables that are missing in CMOR tables
cmore_table['variable_entry']['wo']= OrderedDict([(u'modeling_realm', u'ocean'),
             (u'standard_name', u'sea_water_z_velocity'),
             (u'units', u'm s-1'),
             (u'cell_methods', u'time: mean'),
             (u'cell_measures', u'--OPT'),
             (u'long_name', u'Sea Water Z Velocity'),
             (u'comment',
              u'Not standard CMORE variable'),
             (u'dimensions', u'longitude latitude olevel time'),
             (u'out_name', u'wo'),
             (u'type', u'real'),
             (u'positive', u''),
             (u'valid_min', u''),
             (u'valid_max', u''),
             (u'ok_min_mean_abs', u''),
             (u'ok_max_mean_abs', u'')])

cmore_table['variable_entry']['wpot']= OrderedDict([(u'modeling_realm', u'ocean'),
             (u'standard_name', u'sea_water_z_velocity'),
             (u'units', u'm s-1'),
             (u'cell_methods', u'time: mean'),
             (u'cell_measures', u'--OPT'),
             (u'long_name', u'Vertical Velocity Potential'),
             (u'comment',
              u'Not standard CMORE variable'),
             (u'dimensions', u'longitude latitude olevel time'),
             (u'out_name', u'wpot'),
             (u'type', u'real'),
             (u'positive', u''),
             (u'valid_min', u''),
             (u'valid_max', u''),
             (u'ok_min_mean_abs', u''),
             (u'ok_max_mean_abs', u'')])

cmore_table_ice['variable_entry']['esithick'] = OrderedDict([(u'modeling_realm', u'seaIce'),
             (u'standard_name', u'effective_sea_ice_thickness'),
             (u'units', u'm'),
             (u'cell_methods', u'area: mean where sea_ice time: mean'),
             (u'cell_measures', u'area: areacella'),
             (u'long_name', u'Effective Sea-ice thickness'),
             (u'comment',
              u'Effective thickness of sea ice (volume divided by grid area as was done in CMIP5)'),
             (u'dimensions', u'longitude latitude time'),
             (u'out_name', u'esithick'),
             (u'type', u''),
             (u'positive', u''),
             (u'valid_min', u''),
             (u'valid_max', u''),
             (u'ok_min_mean_abs', u''),
             (u'ok_max_mean_abs', u'')])

#combine ocean and ice variables
cmore_table['variable_entry'].update(cmore_table_ice['variable_entry']) #

# Map FESOM variables to CMOR variables and provide some
# additional information for conversion
vardir = {}
vardir['temp'] = {}
vardir['temp']['dims'] = '3D'
vardir['temp']['cname'] = 'thetao'

vardir['salt'] = {}
vardir['salt']['dims'] = '3D'
vardir['salt']['cname'] = 'so'

vardir['u'] = {}
vardir['u']['dims'] = '3D'
vardir['u']['cname'] = 'uo'
vardir['u']['rotate_with'] = 'v'

vardir['v'] = {}
vardir['v']['dims'] = '3D'
vardir['v']['cname'] = 'vo'
vardir['v']['rotate_with'] = 'u'

vardir['w'] = {}
vardir['w']['dims'] = '3D'
vardir['w']['cname'] = 'wo'

vardir['wpot'] = {}
vardir['wpot']['dims'] = '3D'
vardir['wpot']['cname'] = 'wpot'

vardir['ssh'] = {}
vardir['ssh']['dims'] = '2D'
vardir['ssh']['cname'] = 'zos'
vardir['ssh']['realm'] = 'ocean'

vardir['area'] = {}
vardir['area']['dims'] = '2D'
vardir['area']['cname'] = 'siconc'
vardir['area']['realm'] = 'seaice'

vardir['hice'] = {}
vardir['hice']['dims'] = '2D'
vardir['hice']['cname'] = 'esithick'
vardir['hice']['realm'] = 'seaice'

vardir['uice'] = {}
vardir['uice']['dims'] = '2D'
vardir['uice']['cname'] = 'siu'
vardir['uice']['realm'] = 'seaice'
vardir['uice']['rotate_with'] = 'vice'

vardir['vice'] = {}
vardir['vice']['dims'] = '2D'
vardir['vice']['cname'] = 'siv'
vardir['vice']['realm'] = 'seaice'
vardir['vice']['rotate_with'] = 'uice'

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

def progressbar(progress_total, progress_passed, year, variable, \
                timestep, level, time):
    
    formated_time = num2date(time[timestep], time.units).strftime('%Y-%m-%d')

    sys.stdout.write('{}\n'.format('Variable: '+variable+\
                                   ', Timestep: '+formated_time+\
                                   ', Level: '+str(level)))

    tdif    = progress_total
    tpassed = progress_passed
    ratio   = tpassed/float(tdif)
    filled  = '=' * int( ratio * 50)
    rest    = '-' * ( 50 - int( ratio * 50) )
    sys.stdout.write('|' + filled+'>'+rest+ '| {:.2f}%'.format(ratio*100))
    sys.stdout.write('\r\033[1A')
    sys.stdout.flush()


# Calculate distances and indeces that will be used for interpolation
distances, inds = pf.create_indexes_and_distances(mesh, lons, lats,\
                                                    k=k, n_jobs=8)

# The main loop
def convertit(year):   
    ifile = os.path.join(path_to_data, ifile_template.format(str(year)))
    ofile = os.path.join(path_to_output, ofile_template.format(str(year)))
    
    print('Open {}'.format(ifile))
    fl = Dataset(ifile)
    fw = Dataset(ofile, mode='w',data_model='NETCDF4_CLASSIC', )

    var2d = 0
    var3d = 0
    for varname in out_vars:
        if vardir[varname]['dims'] == '2D':
            var2d += 1
        elif vardir[varname]['dims'] == '3D':
            var3d += 1
    var3d = var3d*len(levels)*fl.variables['time'].shape[0]
    var2d = var2d*fl.variables['time'].shape[0]
    progress_total  = var3d+var2d 
    progress_passed = 0


    # create dimensions
    fw.createDimension('latitude', lons.shape[0])
    fw.createDimension('longitude', lats.shape[1])
    fw.createDimension('time', None)
    fw.createDimension('depth_coord',  levels.shape[0] )

    lat = fw.createVariable('latitude', 'd', ('latitude'))
    lat.setncatts(noempty_dict(cmore_table['axis_entry']['latitude']))
    lat[:] = lats[:,0].flatten()


    lon = fw.createVariable('longitude', 'd', ('longitude'))
    lon.setncatts(noempty_dict(cmore_table['axis_entry']['longitude']))
    lon[:] = lons[0,:].flatten()

    depth = fw.createVariable('depth_coord','d',('depth_coord'))
    depth.setncatts(noempty_dict(cmore_table['axis_entry']['depth_coord']))
    depth[:] = levels

    time = fw.createVariable('time','d',('time'))
    time.setncatts(cmore_table['axis_entry']['time'])
    if distribute_timesteps:
        nsteps = fl.variables['time'].shape[0] 
        td = datetime.timedelta(days = 365/nsteps)
        sdate = datetime.datetime(year,1,1,0,0,0)
        seconds = []
        for i in range(1,nsteps+1):
            workdate = sdate + td*i
            seconds.append( (workdate-sdate).total_seconds() )
        time.units = 'seconds since {}-01-01 00:00:00'.format(str(year))
        time[:] = seconds

    elif fl.variables['time'].units.strip().startswith('seconds since'):
        time.units = fl.variables['time'].units
        time[:] = fl.variables['time'][:]
    elif fl.variables['time'].shape[0] == 12:
        sdate = datetime.datetime(year,1,1,0,0,0)
        td = datetime.timedelta(days = 14.5)
        seconds = []
        for i in range(1,13):
            workdate = datetime.datetime(year,i,1,0,0,0)+td
            seconds.append( (workdate-sdate).total_seconds() )
        time.units = 'seconds since {}-01-01 00:00:00'.format(str(year))
        time[:] = seconds
    else:
        time.units = 'seconds since {}-01-01 00:00:00'.format(str(year))
        time[:] = fl.variables['time'][:]
    

    # Store processed variables (to not repeat 
    # processing for vector variables)
    completed = [] 
    # variables loop
    for varname in out_vars:
        
        # check if we have to convert two variables at once
        # for vector variables.
        do_two_vars = (vardir[varname].has_key('rotate_with') is True) 
                
        #print("Converting {}.".format(varname))
        
        # skip if the variable was already converted
        if varname in completed:
            pass
        # 3D variables processing
        elif vardir[varname]['dims']=='3D':
            
            # Create netCDF variable
            temp = fw.createVariable(vardir[varname]['cname'],'d',\
                                     ('time','depth_coord','latitude','longitude'), \
                                     fill_value=-99999, zlib=zlib, complevel=1)
            # add CMOR complient attributes
            temp.setncatts(noempty_dict(cmore_table['variable_entry'][vardir[varname]['cname']]))
            
            # If we have two convert two variables at once, create netCDF variable for 
            # the second variable
            if do_two_vars is True:
                
                varname2 = vardir[varname]['rotate_with']
                temp2 = fw.createVariable(vardir[varname2]['cname'],'d',('time','depth_coord','latitude','longitude'), fill_value=-99999, zlib=zlib, complevel=1)
                temp2.setncatts(noempty_dict(cmore_table['variable_entry'][vardir[varname2]['cname']]))
                
            # Loop over timesteps for 3D variables
            for i in range(fl.variables[varname].shape[0]):
            #for i in range(2):
                # Get the whole 3D field in to memory. It turns out that this is more 
                # effective than to select individual levels from the file located on the disk. 
                all_layers = fl.variables[varname][i,:]

                # Get the data for the second variable if needed
                if do_two_vars is True:
                    #print("Also converting {}, triggered by {}.".format(varname2, varname))
                    all_layers2 = fl.variables[varname2][i,:]
                # Loop over vertical levels    
                for dlev, llev in enumerate(levels):
                    # get indeces of the gridpoints that corespond to the level
                    ind_depth, ind_noempty, ind_empty = pf.ind_for_depth(llev, mesh)

                    # get the data for the level
                    level_data=np.zeros(shape=(mesh.n2d))
                    level_data[ind_noempty]=all_layers[ind_depth[ind_noempty]]
                    level_data[ind_empty] = np.nan
                    
                    # Spetial treatment of the vector variables that need rotation
                    if do_two_vars is True:

                        # get the data for the level of the second variable
                        level_data2=np.zeros(shape=(mesh.n2d))
                        level_data2[ind_noempty]=all_layers2[ind_depth[ind_noempty]]
                        level_data2[ind_empty] = np.nan
                        
                        #print('Rotate {} and {}'.format(varname, varname2))
                        
                        # Rotate vector variables to geographical grid
                        uunr,vunr = pf.vec_rotate_r2g(angles_for_rotation[0],angles_for_rotation[1], \
                                                angles_for_rotation[2], mesh.x2, mesh.y2,\
                                                level_data, level_data2, 1)
                        
                        # Interpolate rotated variables 
                        #print('interpolation, layer {}'.format(str(llev)))
                        air_nearest = pf.fesom2regular(uunr, mesh, lons, lats, distances=distances,\
                                           inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
                            
                        air_nearest2 = pf.fesom2regular(vunr, mesh, lons, lats, distances=distances,\
                                           inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
                        
                        # Put values to the netCDF variables
                        temp[i,dlev,:,:] = air_nearest[:,:].filled(-99999)
                        temp2[i,dlev,:,:] = air_nearest2[:,:].filled(-99999)
                        
                    else:
                        # Interpolate scalar variable
                        #print('interpolation, layer {}'.format(str(llev)))
                        air_nearest = pf.fesom2regular(level_data, mesh, lons, lats, distances=distances,\
                                           inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
                        # Put values to the netCDF variable
                        temp[i,dlev,:,:] = air_nearest[:,:].filled(-99999)
                    
                    progress_passed += 1
                    
                    if do_two_vars is True:
                        progress_passed += 1
                        
                    progressbar(progress_total, progress_passed, year,\
                                varname, i, llev, time)

            # END Loop over timesteps for 3D variables
            # add variable to the list of processed variables
            completed.append(varname)
            
            if do_two_vars is True:
                completed.append(varname2)
        # End 3D variables processing


        # 2D variables processing
        elif vardir[varname]['dims']=='2D':
            # Create netCDF variable
            temp = fw.createVariable(vardir[varname]['cname'],'d',\
                                     ('time','latitude','longitude'), \
                                     fill_value=-99999, zlib=zlib, complevel=1)

            # add CMOR complient attributes
            temp.setncatts(noempty_dict(cmore_table['variable_entry'][vardir[varname]['cname']]))

            # If we have two convert two variables at once, create netCDF variable for 
            # the second variable
            if do_two_vars is True:
            
                varname2 = vardir[varname]['rotate_with']
                temp2 = fw.createVariable(vardir[varname2]['cname'],'d',\
                                          ('time','latitude','longitude'), \
                                          fill_value=-99999, zlib=zlib, complevel=1)
                temp2.setncatts(noempty_dict(cmore_table['variable_entry'][vardir[varname2]['cname']]))
            
            # For sea ice variables we have to open different file, so
            # open ether ocean or sea ice input file.
            if vardir[varname]['realm']=='ocean':
                temp.setncatts(noempty_dict(cmore_table['variable_entry'][vardir[varname]['cname']]))
                
                ncfile_handler = fl
            
            elif vardir[varname]['realm']=='seaice':
                temp.setncatts(noempty_dict(cmore_table['variable_entry'][vardir[varname]['cname']]))
                
                ifile_ice = os.path.join(path_to_data, ifile_template_ice.format(str(year)))
                ncfile_handler = Dataset(ifile_ice)
                
            # Loop over timesteps for 2D variables  
            for i in range(ncfile_handler.variables[varname].shape[0]):
            #for i in range(2):
                # Get the whole 3D field in to memory. It turns out that this is more 
                # effective than to select individual levels from the file located on the disk. 
                all_layers = ncfile_handler.variables[varname][i,:]
                
                # Get the data for the second variable if needed
                if do_two_vars is True:
                    print("Also converting {}, triggered by {}.".format(varname2, varname))
                    all_layers2 = ncfile_handler.variables[varname2][i,:]
                
                # get indeces of the gridpoints that corespond to the surface level
                ind_depth, ind_noempty, ind_empty = pf.ind_for_depth(0, mesh)

                # get the data for the surface level
                level_data=np.zeros(shape=(mesh.n2d))
                level_data[ind_noempty]=all_layers[ind_depth[ind_noempty]]
                level_data[ind_empty] = np.nan
                
                # Spetial treatment of the vector variables that need rotation
                if do_two_vars is True:
                    # get the data for the surface level of the second variable
                    level_data2=np.zeros(shape=(mesh.n2d))
                    level_data2[ind_noempty]=all_layers2[ind_depth[ind_noempty]]
                    level_data2[ind_empty] = np.nan
                    
                    # Rotate vector variables to geographical grid
                    print('Rotate {} and {}'.format(varname, varname2))
                    uunr,vunr = pf.vec_rotate_r2g(angles_for_rotation[0],angles_for_rotation[1], \
                                            angles_for_rotation[2], mesh.x2, mesh.y2,\
                                            level_data, level_data2, 1)
                    
                    
                    # Interpolate rotated variables )
                    air_nearest = pf.fesom2regular(uunr, mesh, lons, lats, distances=distances,\
                                       inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
                
                    air_nearest2 = pf.fesom2regular(vunr, mesh, lons, lats, distances=distances,\
                                       inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)

                    # fill in netCDF variables
                    temp[i,:,:] = air_nearest[:,:].filled(-99999)
                    temp2[i,:,:] = air_nearest2[:,:].filled(-99999)

                else:
                    # Interpolate scalar variable and fill in netCDF variable.
                    air_nearest = pf.fesom2regular(level_data, mesh, lons, lats, distances=distances,\
                                       inds=inds, radius_of_influence=radius_of_influence, n_jobs=1)
                    temp[i,:,:] = air_nearest[:,:].filled(-99999)
                
                progress_passed += 1

                if do_two_vars is True:
                    progress_passed += 1
               
                progressbar(progress_total, progress_passed, year,\
                                varname, i, 0, time)

            # END Loop over timesteps for 2D variables  
                        
            completed.append(varname)
            
            if do_two_vars is True:
                completed.append(varname2)  


    # end variables loop             

    fw.close()
    print('The {} is ready'.format(ofile))
# end of the main loop
print('='*50)

yearss = range(start_year, end_year+1)
pool = mp.Pool(processes=1)
r = pool.map(convertit, yearss)
pool.close()

#for year in range(start_year, end_year+1):
#    convertit(year)
    # Open input and output netCDF files



