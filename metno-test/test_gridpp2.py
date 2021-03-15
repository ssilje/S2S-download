import xarray as xr
import os
import pandas as pd

import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
import gridpp
import json 
import os

config = {
    'S2S_DIR': '/nird/projects/NS9001K/share/S2S/metno/',
}


fcyear=2020
fcmonth=5
fcday=4

dates_fcycle=pd.date_range(start=f'{fcyear}-{fcmonth}-{fcday}', periods=2, freq='7D') # forecasts start Monday
var_name_abbr='sst'
mdl_vrsn='CY46R1'
S2S_dirbase=config['S2S_DIR']
product='forecast'
curr_date=dates_fcycle[0].strftime('%Y-%m-%d')


def read_grib_file(
    S2S_dirbase,
    product,
    model_version,
    var_name_abbr,
    cast_type,
    date_str,
):
    file_name =  '_'.join([var_name_abbr, model_version, date_str, cast_type, product]) + '.grb'
    file_path = os.path.join(S2S_dirbase, file_name)
    
    print('reading file: ' + file_path)    
    dataopen = xr.open_dataset(file_path, engine='cfgrib')

    return dataopen
  
def make_grid(lats, lons):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_grid = gridpp.Grid(latlats, lonlons)
    return ECMWF_grid 
  
  
  




SAL_GRID1deg = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr='sal',
    cast_type='cf',
    date_str=curr_date,
)

SST_GRID1_5deg = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=curr_date,
)


SST_GRID1_5deg_EUR = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product='forecast_EUR',
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=curr_date,
)



ECMWF_grid1_5deg_EUR = make_grid(SST_GRID1_5deg_EUR.latitude.data, SST_GRID1_5deg_EUR.longitude.data)
ECMWF_grid1_5deg = make_grid(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data)
ECMWF_grid1deg = make_grid(SAL_GRID1deg.latitude.data, SAL_GRID1deg.longitude.data)
