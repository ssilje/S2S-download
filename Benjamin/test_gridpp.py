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

lead_time=np.arange(1,47)
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


ECMWF_grid1_5deg = make_grid(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data)
ECMWF_grid1deg = make_grid(SAL_GRID1deg.latitude.data, SAL_GRID1deg.longitude.data)


with open("metadata_BW_sites.json") as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
BW_grid = gridpp.Points(data_BW.lat, data_BW.lon)


step = SST_GRID1_5deg.get_index('step')[0]
#for step in ds_cf.get_index('step'):
print('Time step: ' + str(step))

### Interpolating to 1 deg to make it work
SST_grid1deg = gridpp.bilinear(
   ECMWF_grid1_5deg,
   ECMWF_grid1deg,
   np.transpose(SST_GRID1_5deg.sst[step.days - 1,:,:].data)
)
     
SST_1deg2point = gridpp.bilinear(
    ECMWF_grid1deg, 
    BW_grid, 
    gridpp.fill_missing(SST_grid1deg)
)

SST_1_5deg2point = gridpp.bilinear(
    ECMWF_grid1_5deg, 
    BW_grid, 
    gridpp.fill_missing(np.transpose(SST_GRID1_5deg.sst[step.days - 1,:,:].data))
)
    
