import click
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Qt4Agg')
#%matplotlib inline
import matplotlib.pylab as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cmocean import cm as cmo
from matplotlib import cm
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
import pyfesom as pf
from cartopy.util import add_cyclic_point


@click.command()
@click.argument('meshpath', type=click.Path(exists=True))
@click.argument('ifile', type=click.Path(exists=True))
@click.argument('variable', default='temp', required=False)
@click.option('--depth', '-d', default=0, type=click.FLOAT, show_default=True,
               help='Depth in meters.')
@click.option('--box', '-b',
              nargs=4,
              type=(click.IntRange(-180, 180),
                    click.IntRange(-180, 180),
                    click.IntRange(-90, 90),
                    click.IntRange(-90, 90)),
              default=(13, 30, 54, 66), show_default=True,
              help='Map boundaries in -180 180 -90 90 format.')
@click.option('--timestep', '-t', default=0, show_default=True,
              help='Timstep from netCDF variable, strats with 0.')
@click.option('--minmax', '-l', nargs=2, type=click.INT,
              help='Minimun and maximum values for plotting.')
@click.option('--mapproj','-m', type=click.Choice(['merc', 'pc','np','sp', 'rob']), default='merc', help = 'Map projection. Options are Mercator (merc), Plate Carree (pc), North Polar Stereo (np), South Polar Stereo (sp),  Robinson (rob)')
@click.option('--quiet', '-q', is_flag=True,
              help='If present additional information will not be printed.')
@click.option('--ofile','-o' ,type=click.Path(exists=False) ,
              help='Path to the output figure. If present the image will be saved to the file instead of showing it. ')
def showo(meshpath, ifile, variable, 
          depth, box, timestep, minmax, 
          mapproj, quiet, ofile):
    '''
    meshpath - Path to the folder with FESOM1.4 mesh files.

    ifile    - Path to FESOM1.4 netCDF file.

    variable - The netCDF variable to be plotted.
    '''
    if not quiet:
        click.secho('Mesh: {}'.format(meshpath))
        click.secho('File: {}'.format(ifile))
        click.secho('Variable: {}'.format(variable), fg='red')
        click.secho('Depth: {}'.format(depth), fg='red')
        click.secho('BOX: {}'.format(box))
        click.secho('Timestep: {}'.format(timestep))
        if minmax:
            click.secho('Min/max: {}'.format(minmax), fg='red')
        else:
            click.secho('Min/max: auto', fg='red')

    mesh = pf.load_mesh(meshpath)
    flf = Dataset(ifile)
    left, right, down, up = box

    elem_nonan, no_nan_tri = pf.cut_region(mesh, [left, right, down, up], 0)
    if variable == 'topo':
        level_data = mesh.topo
    else:
        level_data , nnn = pf.get_data(flf.variables[variable][timestep,:], mesh, depth)
    
    if minmax:
        mmin, mmax = minmax
    else:
        mmin = level_data[elem_nonan].min()
        mmax = level_data[elem_nonan].max()
    
    plt.figure(figsize=(10,10))
    if mapproj == 'merc':
        ax = plt.subplot(111, projection=ccrs.Mercator())
    elif mapproj == 'pc':
        ax = plt.subplot(111, projection=ccrs.PlateCarree())
    elif mapproj == 'np':
        ax = plt.subplot(111, projection=ccrs.NorthPolarStereo())
    elif mapproj == 'sp':
        ax = plt.subplot(111, projection=ccrs.SouthPolarStereo())
    elif mapproj == 'rob':
            ax = plt.subplot(111, projection=ccrs.Robinson())

    ax.set_extent(box, crs=ccrs.PlateCarree())
    mm = ax.tripcolor(mesh.x2, 
                      mesh.y2,
                      elem_nonan, 
                      level_data,
                      transform=ccrs.PlateCarree(), 
                      edgecolors='k',
                      cmap=cm.Spectral_r,
                      vmin=mmin,
                      vmax=mmax
                      )
    ax.coastlines(lw=0.5, resolution='10m')
    plt.colorbar(mm, orientation='horizontal', pad =0.03, )
    plt.title('{} at {}'.format(variable, depth), size = 20)
    if ofile:
        plt.savefig(ofile, dpi=200)
    else:
        plt.show()
if __name__ == '__main__':
    showo()
