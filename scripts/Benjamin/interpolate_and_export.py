import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import gridpp
import json 
import os

from S2S.file_handling import read_grib_file
from S2S.local_configuration import config
from S2S.gridpp_helpers import make_grid_object

# Dates
# var_name='sav300' 

lead_time=np.arange(1,47)
fc_year=2020
fc_month=5
fc_day=4

dates_fcycle=pd.date_range(start=f'{fc_year}-{fc_month}-{fc_day}', periods=2, freq='7D') # forecasts start Monday

# Read in data for a given date
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

# Grid/points for interpolation from/to

ECMWF_grid = make_grid_object(ds_cf)
ECMWF_grid1deg = make_grid_object(grid1deg)

with open(os.path.join(config['BW_DIR'], "metadata_BW_sites.json")) as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
BW_grid = gridpp.Points(data_BW.lat, data_BW.lon)

#
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
