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


## funkar ikkje
def make_grid_flatten(lats, lons):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_grid = gridpp.Grid(latlats.flatten(), lonlons.flatten())
    return ECMWF_grid 
  
  
def make_points(lats, lons,valid_points):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    print(latlats.shape)
    print(lonlons.shape)
    print(valid_points.shape)
    #valid_points=np.transpose(valid_points)
    print(valid_points.shape)
    print(latlats[valid_points].shape)

#    ECMWF_points = gridpp.Points(latlats[valid_points], lonlons[valid_points])
    ECMWF_points = gridpp.Points(latlats, lonlons)    
    return ECMWF_points   

def make_points_flatten(lats, lons,valid_points):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_points = gridpp.Points(latlats.flatten()[valid_points], lonlons.flatten()[valid_points])
    return ECMWF_points


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



SST_GRID1_5deg_EUR = SST_GRID1_5deg 


with open("metadata_BW_sites.json") as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
BW_grid = gridpp.Points(data_BW.lat, data_BW.lon)
print(data_BW.lat.shape)
ECMWF_grid1_5deg_EUR = make_grid(SST_GRID1_5deg_EUR.latitude.data, SST_GRID1_5deg_EUR.longitude.data)
ECMWF_grid1_5deg = make_grid(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data)
ECMWF_grid1deg = make_grid(SAL_GRID1deg.latitude.data, SAL_GRID1deg.longitude.data)
#ECMWF_grid1_5deg_EUR = make_points(SST_GRID1_5deg_EUR.latitude.data, SST_GRID1_5deg_EUR.longitude.data)
#ECMWF_grid1_5deg_EUR = make_grid_flatten(SST_GRID1_5deg_EUR.latitude.data, SST_GRID1_5deg_EUR.longitude.data)
SST_1_5deg2point = gridpp.nearest(
    ECMWF_grid1_5deg, 
    BW_grid, 
    gridpp.fill_missing(np.transpose(SST_GRID1_5deg.sst[0,:,:].data))
)
print(SST_1_5deg2point)

sst = SST_GRID1_5deg_EUR.sst[0,:,:].data.flatten()
valid_points = np.isnan(sst) == 0  # Ocean points

ECMWF_points_EUR = make_points_flatten(SST_GRID1_5deg_EUR.latitude.data, SST_GRID1_5deg_EUR.longitude.data, valid_points) 


ECMWF_points_EUR_v2 = gridpp.Points(ECMWF_grid1_5deg_EUR.get_lats().flatten()[valid_points],
     ECMWF_grid1_5deg_EUR.get_lons().flatten()[valid_points])

#SST_points_EUR_v2 = gridpp.nearest(ECMWF_points_EUR_v2, BW_grid, sst[valid_points])


#SST_points_EUR = gridpp.nearest(ECMWF_points_EUR, BW_grid, sst[valid_points])



sst = SST_GRID1_5deg_EUR.sst[0,:,:].data
sst=np.transpose(sst)
valid_points = np.isnan(sst) == 0  # Ocean points



ECMWF_points_EUR = make_points(SST_GRID1_5deg_EUR.latitude.data, SST_GRID1_5deg_EUR.longitude.data, valid_points)
# valid_points=np.transpose(valid_points)
 
#SST_points_EUR = gridpp.nearest(ECMWF_points_EUR, BW_grid, sst[valid_points])        
print(valid_points.shape)
print(sst.shape)
print((sst[valid_points]).shape)
SST_points_EUR = gridpp.nearest(ECMWF_points_EUR, BW_grid, sst) 
