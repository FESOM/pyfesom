import glob
from joblib import Parallel, delayed
import os
import click


def cdo_command(ifile, opath, command, ext, skip):
    if opath != 'no':
        ofile = os.path.join(opath, '{}_{}.nc'.format(os.path.basename(ifile)[:-3], ext))
    else:
        ofile = ' '

    if skip:
        if os.path.isfile(ofile):
            print('File {} exist, --checko flag is present, skipping'.format(ofile))
            return

    print('cdo {} {} {}'.format(command, ifile, ofile))
    os.system('cdo {} {} {}'.format(command, ifile, ofile))

@click.command()
@click.argument('ipath', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('opath', nargs=1, required=False, default='no')
@click.option('--ncore', '-n', default=2, help = 'Number of cores (parallel processes)', show_default=True)
@click.option('--cdo','-c', required=True, help = 'CDO command as a string !!!IN QUOTATION MARKS!!!, eg \" monmean -shifttime,-12hour \"')
@click.option('--ext','-e', default='tm', required=False,show_default=True,
              help='Extention to be used for the output file.')
@click.option('--skip', '-s', is_flag=True,
              help='Skip the calculation if the output file already exist.')
def pcdo(ipath, opath, ncore, cdo, ext, skip):
    '''
    Runs several (-n) cdo processes in paralel. Input (ipath) is a list (wildcard) of files. The cdo command (-c) is
    executed for every file and the output files with extention (-e) will be written to the output path (opath).

    Example:

    python pcdo.py

    ipath - Input files, must be the path with wildcards (e.g. /path/to/files/temp_fesom_193[3-7]????.nc)

    opath - Path where the output will be stored or "no" for operators that do not require output file.
    '''
    Parallel(n_jobs=ncore)(delayed(cdo_command)(l, opath, cdo, ext, skip) for l in ipath)

if __name__ == '__main__':
    pcdo()


