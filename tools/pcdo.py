import glob
from joblib import Parallel, delayed
import os
import click


def cdo_comand(ifile, opath, command, ext, checko):
    if opath != 'no':
        ofile = os.path.join(opath, '{}_tm.nc'.format(os.path.basename(ifile)[:-3]))
    else:
        ofile = ' '

    if checko:
        if os.path.isfile(ofile):
            print('File {} exist, --checko flag is present, skipping'.format(ofile))
            return

    print('cdo {} {} {}'.format(command, ifile, ofile))
    os.system('cdo {} {} {}'.format(command, ifile, ofile))

@click.command()
@click.argument('ipath', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('opath', nargs=1, required=False, default='no')
@click.option('--ncore', '-n', default=2, help = 'Number of cores (parallel processes)', show_default=True)
@click.option('--cdo', required=True, help = 'String of cdo commands !!!IN QUOTATION MARKS!!!, eg \" monmean -shifttime,-12hour \"')
@click.option('--ext','-e', default='tm', required=False,show_default=True,
              help='Extention to be used for the output file.')
@click.option('--checko', '-c', is_flag=True, 
              help='Skip the calculation if the output file already exist.')
def pcdo(ipath, opath, ncore, cdo, ext, checko):
    '''
    ipath - Input files, must be the path with wildcards (e.g. /path/to/files/temp_fesom_193[3-7]????.nc)

    opath - Path where the output will be stored or "no" for operators that do not require output file.
    '''
    Parallel(n_jobs=ncore)(delayed(cdo_comand)(l, opath, cdo, ext, checko) for l in ipath)

if __name__ == '__main__':
    pcdo()


