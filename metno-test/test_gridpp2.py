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
  

def make_points_flatten(lats, lons,valid_points):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_points = gridpp.Points(latlats.flatten()[valid_points], lonlons.flatten()[valid_points])
    return ECMWF_points

def make_grid_flatten(lats, lons,valid_points):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_points = gridpp.Grid(latlats.flatten()[valid_points], lonlons.flatten()[valid_points])
    return ECMWF_points

SST_GRID1_5deg = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=curr_date,
)



with open("metadata_BW_sites.json") as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
    
BW_grid = gridpp.Points(data_BW.lat, data_BW.lon)

ECMWF_grid1_5deg = make_grid(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data)

SST_1_5deg2point = gridpp.nearest(
    ECMWF_grid1_5deg, 
    BW_grid, 
    gridpp.fill_missing(np.transpose(SST_GRID1_5deg.sst[0,:,:].data))
)
print('test 1 med grid2points')
print(SST_1_5deg2point)

sst = SST_GRID1_5deg.sst[0,:,:].data.flatten()
sst = np.transpose(sst)
valid_points = np.isnan(sst) == 0  # Ocean points

ECMWF_points = make_points_flatten(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data, valid_points) 
 
SST_points = gridpp.nearest(ECMWF_points, BW_grid, sst[valid_points])        
print('test 2 med points2points')
print(SST_points)


#ECMWF_grid1_5deg_valid_points = make_grid_flatten(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data, valid_points)
#SST_grid_points = gridpp.nearest(ECMWF_grid1_5deg_valid_points, BW_grid, sst[valid_points])
