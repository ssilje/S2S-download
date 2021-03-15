#%%
import pandas as pd
import numpy as np
import gridpp
import json
import os

from load_functions import read_grib_file, make_grid
from local_configuration import config
from interpolate_and_export import interpolate_and_save

#%% Fill in instance data

var_name='sst'
var_name_abbr='sst'
mdl_vrsn='CY46R1'
S2S_dirbase=config['S2S_DIR']
product='forecast'

# sst
lead_time=np.arange(1,47)
fcyear=2019
fcmonth=7
fcday=1

# sav300
# fcyear=2020
# fcmonth=5
# fcday=4


dates_fcycle=pd.date_range(start=f'{fcyear}-{fcmonth}-{fcday}', periods=2, freq='7D') # forecasts start Monday
date_str=dates_fcycle[0].strftime('%Y-%m-%d')

#%% Read in data for a given date

ds_cf = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=date_str,
)

ds_pf = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='pf',
    date_str=date_str,
)

#%%
ECMWF_grid = make_grid(ds_cf.latitude.data, ds_cf.longitude.data)

with open(os.path.join(config['BW_DIR'], "metadata_BW_sites.json")) as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
BW_points = gridpp.Points(data_BW.lat, data_BW.lon)

#%%
interpolate_and_save(
    var_name=var_name,
    product=product,
    date_str=date_str,
    ds_cf=ds_cf,
    ds_pf=ds_pf,
    in_grid=ECMWF_grid,
    out_points=BW_points,
    data_export_dir=config['CF_DATA_DIR'],
    data_BW=data_BW,
)
#%%

# %%
