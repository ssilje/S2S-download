import xarray as xr
import numpy as np
import pandas as pd
import json
import netCDF4
import h5py

from datetime import datetime,timedelta

# TODO: get setup from config
#       append to exisiting local file

# local dependencies
from .nird2local import pull, read
from .grid import cap_grid
from .handle_domain import update_storage_dict, get_bounds

# setup
username   = 'heau'
uni_path   = '/SFE/Systems_monthly_grib'
var_name   = '/sea_surface_temperature'
file_c     = '.grib'
domainID   = 'NVK'
ftype      = 'hindcast'
model      = 'ecmwf'
var_abb    = 'sst'
t_res      = 'monthly'
i_time     = (1993,1) # start time, fmt: (year,month)
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
step       = np.array([])
number     = np.array([])
out_var    = np.array([])
steps      = np.arange(12*(e_time[0]-i_time[0]))

# generate list of filenames to read
filenames = [
                '_'.join([var_name,model,\
                str(i_time[0]+ii//12),\
                str(i_time[1]+ii%12)])\
                + file_c for ii in steps
            ]

# copy data from nird
local_path = pull(uni_path,filenames,username)

# get domain boundaries
bounds = get_bounds(domainID)

# downsize grid - store locally
for n,filename in enumerate(filenames):

    data = xr.open_dataset(local_path+filename, engine='cfgrib')

    lon,lat,var = cap_grid(
                            lon=data['longitude'],
                            lat=data['latitude'],
                            var=data[var_abb],
                            bounds=bounds
                            )
    step   = data['step']
    number = data['number']
    time   = data['time']


    print(filename)
print(data)
