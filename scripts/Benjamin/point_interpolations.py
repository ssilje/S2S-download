#%%
import pandas as pd
import numpy as np
import xarray as xr
import xarray_extras as xr_e
import gridpp
import json
import os

from S2S.date_helpers import get_forcast_date_cycle
from S2S.gridpp_helpers import make_grid_from_grb, make_points_from_grb
from S2S.file_handling import read_grib_file
from S2S.local_configuration import config
#%%
dates_fcycle = get_forcast_date_cycle(
    start_year=2020,
    start_month=5,
    start_day=4,
    num_weeks=1,
)

#%%
mdl_vrsn = 'CY46R1'
dirbase = config['S2S_DIR']
product = 'forecast'
curr_date = dates_fcycle[0]


with open(os.path.join(config['BW_DIR'], 'sites.json')) as json_file:
    data_BW = pd.DataFrame(json.load(json_file))

out_points = gridpp.Points(data_BW.lat, data_BW.lon)
out_IDs = data_BW.localityNo

#%% Define objects that are reused
grb_data_sst_dummy = read_grib_file( dirbase=dirbase, product=product, model_version=mdl_vrsn, var_name_abbr='sst', cast_type='cf', date=dates_fcycle[0])
valid_points_sst = np.transpose(np.isnan(grb_data_sst_dummy.variables['sst'].isel(step = 0).data) == 0) # NB: Gridpp and xarray have opposite lat-lon axes
in_points_sst = make_points_from_grb(grb_data_sst_dummy, valid_points_sst)

grb_data_sav300_dummy = read_grib_file( dirbase=dirbase, product=product, model_version=mdl_vrsn, var_name_abbr='sal', cast_type='cf', date=dates_fcycle[0])
valid_points_sav300 = np.transpose(np.isnan(grb_data_sav300_dummy.variables['sav300'].isel(step = 0).data) == 0) # NB: Gridpp and xarray have opposite lat-lon axes
in_points_sav300 = make_points_from_grb(grb_data_sav300_dummy, valid_points_sav300)

# Initialize empty xarray to insert in
data_cf_empty = xr.DataArray(dims = ('step', 'locNo', 'date'), 
    coords={
        'step': grb_data_sst_dummy.get_index('step'),
        'locNo': out_IDs,
        'date': dates_fcycle,
    }
)
data_converted_cf = xr.Dataset(
    {
        'sst': data_cf_empty.copy(), 
        'sav300': data_cf_empty.copy(), 
    }
)

# Perform interpolations for each date and step
for curr_date in dates_fcycle:
    # Sea surface temperature
    grb_data_sst = read_grib_file(dirbase=dirbase, product=product, model_version=mdl_vrsn, var_name_abbr='sst', cast_type='cf', date=curr_date)
    for curr_step in grb_data_sst_dummy.get_index('step'):
        in_values_sst = np.transpose(grb_data_sst.sel(step = curr_step).variables['sst'].data)[valid_points_sst]
        out_values_sst = gridpp.nearest(
            in_points_sst,
            out_points,
            in_values_sst
        )
        data_converted_cf.sel(step = curr_step, date = curr_date).variables['sst'].data[:] = out_values_sst

    # Salinity
    grb_data_sav300 = read_grib_file(dirbase=dirbase, product=product, model_version=mdl_vrsn, var_name_abbr='sal', cast_type='cf', date=curr_date)
    for curr_step in grb_data_sst_dummy.get_index('step'):
        in_values_sav300 = np.transpose(grb_data_sav300.sel(step = curr_step).variables['sav300'].data)[valid_points_sav300]
        out_values_sav300 = gridpp.nearest(
            in_points_sav300,
            out_points,
            in_values_sav300
        )
        data_converted_cf.sel(step = curr_step, date = curr_date).variables['sav300'].data[:] = out_values_sav300
    
#%% Save to file, with and without stepping for different uses

data_converted_cf.to_dataframe().to_csv(f'data_wsteps_cf_{dates_fcycle[0].date()}_{dates_fcycle[-1].date()}.csv')
data_converted_cf.isel(step = 0).drop('step').to_dataframe().to_csv(f'data_nsteps_cf_{dates_fcycle[0].date()}_{dates_fcycle[-1].date()}.csv')

# %%
