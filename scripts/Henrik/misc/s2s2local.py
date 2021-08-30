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
var_name   = '/sea_surface_temperature'
file_c     = '.grib'
domainID   = 'NVK'
ftype      = 'forecast'
model      = 'ecmwf'
var_abb    = 'sst'
t_res      = 'monthly'
cast_type  = 'cf'
model_v    = 'CY46R1_CY47R1'
i_time     = (2020,1,23) # start time, fmt: (year,month)
e_time     = (2021,1,14) # end time, fmt: (year,month)
# dtype      = 'float64'
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
uni_path   = '/S2S/Benjamin/'+ftype+'/ECMWF/sfc/sst/'


dates_fcycle_1 = pd.date_range(
                        start=f'{i_time[0]}-{i_time[1]}-{i_time[2]}',
                        end=f'{e_time[0]}-{e_time[1]}-{e_time[2]}',
                        freq='7D'
                        )

dates_fcycle_2 = pd.date_range(
                        start=f'{i_time[0]}-{i_time[1]}-{i_time[2]+4}',
                        end=f'{e_time[0]}-{e_time[1]}-{e_time[2]+4}',
                        freq='7D'
                        )

for dates_fcycle in [dates_fcycle_1,dates_fcycle_1]:

    filenames = ['_'.join([
                            var_abb,
                            model_v,
                            date.strftime('%Y-%m-%d'),
                            cast_type,
                            ftype
                            ])
                        +'.grb'
                        for date in dates_fcycle
                        ]

    # copy data from nird
    local_path = pull(uni_path,filenames,username)
