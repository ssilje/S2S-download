import xarray as xr
import numpy as np
import pandas as pd
import json
import netCDF4
import h5py

# TODO: get setup from config
#       append to exisiting local file

# local dependencies
from .nird2local import pull, read
from .grid import cap_grid
from .handle_domain import update_storage_dict, get_bounds

dir_copy = 1

# setup
username   = 'heau'
uni_path   = '/SFE/ERA_monthly_nc'
var_name   = '/sea_surface_temperature'
file_c     = '.nc'
domainID   = 'NVK'
ftype      = 'reanalysis'
model      = 'era'
var_abb    = 'sst'
t_res      = 'monthly'
i_time     = (2020,1) # start time, fmt: (year,month)
e_time     = (2020,9) # end time, fmt: (year,month)
dtype      = 'float64'
store_path = '_'.join(['./data/EIDE/'+model,ftype,var_abb,t_res,domainID,
                str(i_time[0]),str(i_time[1])+'-'+str(e_time[0]),str(e_time[1])])\
                +'.hdf5'

# initialize variables
filenames  = []
filename   = ''
bounds     = ()
data       = None # netCDF
lon        = np.array([])
lat        = np.array([])
var        = np.array([])
out_var    = np.array([])
steps      = np.arange(12*(e_time[0]-i_time[0]))

# generate list of filenames to read
filenames = [
                '_'.join([var_name,model,\
                str(i_time[0]+ii//12),\
                str(i_time[1]+ii%12)])\
                for ii in steps
            ]

if dir_copy:

    bounds = get_bounds(domainID)

    # re-initialize out_var
    data = read(uni_path,filenames[0]+file_c,username)
    lon,lat,var = cap_grid(
                            lon=data['longitude'],
                            lat=data['latitude'],
                            var=data[var_abb],
                            bounds=bounds
                            )
    out_var = np.empty((steps.shape[0],lat.shape[0],lon.shape[0]))

    # read files from nird
    for n,filename in enumerate(filenames):

        data = read(uni_path,filename+file_c,username)

        lon,lat,var = cap_grid(
                                lon=data['longitude'],
                                lat=data['latitude'],
                                var=data[var_abb],
                                bounds=bounds
                                )
        out_var[n,:,:] = var.squeeze()

        print(filename)

    # store files
    with h5py.File(store_path, 'w') as h5py_file:

        h5py_file.create_dataset(var_abb, data=out_var, dtype='float64')
        h5py_file.create_dataset('lon', data=lon, dtype='float64')
        h5py_file.create_dataset('lat', data=lat, dtype='float64')

        h5py_file.close()

    update_storage_dict(
                        domainID,
                        ftype,
                        model,
                        var_abb,
                        t_res,
                        i_time,
                        e_time,
                        dtype,
                        store_path
                        )

else:

    local_path = pull(uni_path,filenames,username)
