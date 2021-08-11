#### rename the high resolution hindcast files that currently have two model versions in the file string to contain only the valid model version string

# load required packages
import subprocess as sp
import glob
import re
from datetime import datetime

# load own functions
exec(open('/nird/projects/NS9001K/owul/projects/s2s_global_skill/scripts/fcts.py').read())

# write output to file (inside current working directory!)
with open("log_rename_modelversion.txt", "a") as o:
    o.write('\n\n\t---execution at {:}---\n'.format(datetime.today()))
    
    # input to identify the correct path:
    MODEL = 'ECMWF'
    hcfc = 'hindcast' # 'hindcast' 'forecast'
    hres = False # True False

    if hres:
        hresext = '/MARS'
        resol = 'high resolution'
    else:
        hresext = ''
        resol = 'low resolution'

    # location of files:
    fc_path = '/nird/projects/NS9853K/DATA/S2S{0:s}/{1:s}/{2:s}/sfc/'.format(hresext,hcfc,MODEL)

    o.write('\nINPUT:\n\t{1:s}\n\t{2:s}\n\t{3:s}\
    \nworking in path {0:s}\n'.format(fc_path,MODEL,hcfc,resol))

    allvar_path = glob.glob(fc_path+'*/')

    # write the output of the following to file:
    for var_path in allvar_path:
        # get all .grb file names from the directory:
        file_list = glob.glob(var_path + '*.grb')

        # loop through all of them
        for fi in file_list:
            # extract the date string from the file name by searching for the date pattern:
            datestring = re.search('\d\d\d\d-\d\d-\d\d',fi.split('/')[-1]).group()

            # find the correct model version matching the datec using function from fcts.py 
            valid_model_version = which_mv_for_init(datestring,MODEL)

            if valid_model_version is not None:
                # different treatment when no model version is part of the filename than if it has one or more:
                if len(fi.split('CY')) >= 2:
                    # define new file name with just the single model version:
                    pattern_MV = ''
                    for ii in range(len(fi.split('CY'))-1):
                        pattern_MV += '_CY\d\dR\d'
                    
                    # get the exact sequence from the file string matching the pattern:
                    mv_long_string = re.search(pattern_MV,fi).group()
                    
                    fi_out = fi.replace(mv_long_string,'_{0:s}'.format(valid_model_version))
                elif len(fi.split('CY')) < 2:
                    # for high resolution files, the model version comes before the resolution specification
                    if 'MARS' in fc_path:
                        fi_out = fi.replace('_05x05_','_{0:s}_05x05_'.format(valid_model_version))
                    # for low resolution files, there is no grid specification and the model version comes right before the date:
                    else:
                        fi_out = fi.replace('_{0:s}_'.format(datestring),'_{1:s}_{0:s}_'.format(datestring,valid_model_version))

                hcfc_extension = 'f_{0:s}.grb'.format(hcfc)
                # if the extension exists, remove it:
                if hcfc_extension in fi: # note: if it folder name and extension, do not match, the extension is kept!
                    fi_out = fi_out.replace(hcfc_extension,'f.grb')

                if (fi != fi_out) and (datestring == re.search('\d\d\d\d-\d\d-\d\d',fi_out.split('/')[-1]).group()):
                    o.write('\nmv -u\n{0:s}\n{1:s}\n'.format(fi,fi_out))
                    ####### execute the actual renaming with a shell command ######
                    # run mv with -u option to keep the most current file in case both exist 
                    sp.call('mv -u {0:s} {1:s}'.format(fi,fi_out),shell=True)
                elif (fi == fi_out):
                    o.write('\nwill not rename {0:s},\nfilename seems to be correct\n'.format(fi))
                elif (datestring != re.search('\d\d\d\d-\d\d-\d\d',fi_out.split('/')[-1]).group()):
                        o.write('\nwill not rename {0:s},\ndates do not match with {1:s}\n'.format(fi,fi_out))
            else:
                o.write('\nskipping file for {0:s}, no matching model version found\n'.format(datestring))
                continue