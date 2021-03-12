#%%
import sys  
sys.path.append('/nird/projects/NS9001K/sso102/Python/S2S-download/lib')  

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
import gridpp
import json 
import os


from load_functions import read_grib_file
from local_configuration import config

#%% Dates
# var_name='sav300' 

lead_time=np.arange(1,47)
fcyear=2020
fcmonth=5
fcday=4

dates_fcycle=pd.date_range(start=f'{fcyear}-{fcmonth}-{fcday}', periods=2, freq='7D') # forecasts start Monday

#%% Read in data for a given date

var_name_abbr='sst'
mdl_vrsn='CY46R1'
S2S_dirbase=config['S2S_DIR']
product='forecast'
curr_date=dates_fcycle[1].strftime('%Y-%m-%d')

grid1deg = read_grib_file(
    S2S_dirbase='/nird/projects/NS9001K/sso102/S2S/DATA/',
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr='sal',
    cast_type='cf',
    date_str=curr_date,
)

ds_cf = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=curr_date,
)
ds_pf = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='pf',
    date_str=curr_date,
)

#%% Grid/points for interpolation from/to
def make_grid(lats, lons):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_grid = gridpp.Grid(latlats, lonlons)
    return ECMWF_grid
ECMWF_grid = make_grid(ds_cf.latitude.data, ds_cf.longitude.data)
ECMWF_grid1deg = make_grid(grid1deg.latitude.data, grid1deg.longitude.data)


with open(os.path.join(config['BW_DIR'], "metadata_BW_sites.json")) as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
BW_grid = gridpp.Points(data_BW.lat, data_BW.lon)

#%%
file_path = os.path.join(config['CF_DATA_DIR'], f'{product}_{curr_date}_bilinear_sst.csv')

print('File save location: ' + file_path)

df_out = pd.DataFrame(
    columns = tuple([
        'locNo', 
        'date',
        'step_fwrd', 
        'sst_ctrl'
    ] + [
        'sst_ptrb_%02d'%(i)
        for i in ds_pf.get_index('number')
    ])
)
for step in ds_cf.get_index('step'):
    print('Time step: ' + str(step))
    # NB: Is this the best way to deal with missings from on-land coordinates? 
    # NB: Axes of lat/lon are reversed between gridpp and xarray.
    cf_ovalues_grid1deg = gridpp.bilinear(
       ECMWF_grid,
       ECMWF_grid1deg,
       np.transpose(ds_cf.sst[step.days - 1,:,:].data)
    )
     
    cf_values = gridpp.bilinear(
        ECMWF_grid1deg, 
        BW_grid, 
        gridpp.fill_missing(cf_ovalues_grid1deg)
    )
    
 #   cf_nearest= gridpp.nearest(
 ##      ECMWF_grid,
  #     BW_grid, 
   #    np.transpose(ds_cf.sst[step.days - 1,:,:].data)
   # )
    

    pf_values = np.empty((len(data_BW), len(ds_pf.get_index('number'))), dtype=float)
    for num in ds_pf.get_index('number'):
        pf_ovalues_grid1deg = gridpp.bilinear(
        ECMWF_grid,
        ECMWF_grid1deg,
        np.transpose(ds_pf.sst[num - 1, step.days - 1,:,:].data)
       )
            
        pf_values[:, num - 1] = gridpp.bilinear(
            ECMWF_grid1deg, 
            BW_grid, 
            gridpp.fill_missing(pf_ovalues_grid1deg)
        )
    for i in range(len(data_BW)):
        row_dict = {
            'locNo' : data_BW['localityNo'][i],
            'step_fwrd' : step.days,
            'date' : curr_date,
            'sst_ctrl' : cf_values[i],
        }
        for num in ds_pf.get_index('number'):
            row_dict['sst_ptrb_%02d'%(num)] = pf_values[i, num - 1]

        df_out = df_out.append(pd.Series(row_dict), ignore_index = True)
df_out.to_csv(file_path)


# %%
