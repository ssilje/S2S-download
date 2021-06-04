
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import time as time_lib

from S2S.local_configuration import config
from S2S.data_handler        import BarentsWatch, ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.graphics.graphics as gr

domainID = 'norwegian_coast'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

process_hindcast     = True
process_observations = True
calibrate            = True
verify               = True

observations = BarentsWatch().load(['Hisdalen','Langoy S']).sortby('time')

print('Process hindcast')
if process_hindcast:

    print('\tLoad hindcast')
    hindcast = ECMWF_S2SH(high_res=False)\
                    .load(var,t_start,t_end,domainID)[var]-272.15

    print('\tInterpolate hindcast to point locations')
    hindcast = hindcast.sortby(['time','step'])\
                    .interpolate_na(
                                    dim='lon',
                                    method='nearest',
                                    fill_value="extrapolate"
                                )\
                        .interp(
                                lon=observations.lon,
                                lat=observations.lat,
                                method='nearest'
                            )
    print('\tApply 7D running mean')
    hindcast = hindcast\
                    .rolling(step=7,center=True).mean()\
                        .dropna('step')

    print('\tAssign validation time')
    hindcast = xh.assign_validation_time(hindcast)

    print(hindcast)

    hindcast.to_netcdf(config['VALID_DB']+'/h_temp.nc')

    process_observations = True

hindcast = xr.open_dataset(config['VALID_DB']+'/h_temp.nc')

print('Process observations')
if process_observations:

    clim_mean,clim_std = xh.o_climatology(observations[var],'time.month')

    val_obs   = xh.at_validation(observations,hindcast.validation_time)
    clim_mean = xh.at_validation(clim_mean,hindcast.validation_time,ddays=3.5)
    clim_std  = xh.at_validation(clim_std,hindcast.validation_time,ddays=3.5)

    val_obs.to_netcdf(config['VALID_DB']+'/v_temp.nc')
    clim_mean.to_netcdf(config['VALID_DB']+'/clim_mean_temp.nc')
    clim_std.to_netcdf(config['VALID_DB']+'/clim_std_temp.nc')

    model  = True

val_obs   = xr.open_dataset(config['VALID_DB']+'/v_temp.nc')
clim_mean = xr.open_dataset(config['VALID_DB']+'/clim_mean_temp.nc')
clim_std  = xr.open_dataset(config['VALID_DB']+'/clim_std_temp.nc')

print('Calibrate')
if calibrate:

    mean,std = xh.c_climatology(hindcast,dim='validation_time.month')

    mean = xh.climatology_to_validation_time(mean,hindcast.time+hindcast.step)
    std  = xh.climatology_to_validation_time(std,hindcast.time+hindcast.step)

    hindcast = (hindcast-mean)/std

    hindcast.to_netcdf(config['VALID_DB']+'/hp_temp.nc')
    verify = True

hindcast = xr.open_dataset(config['VALID_DB']+'/hp_temp.nc')

print('Verify')
if verify:

    hindcast   = (hindcast*clim_std) + clim_mean
    clim_fc    = models.clim_fc(clim_mean,clim_std)
    hindcast   = hindcast.sel(time=slice(val_obs.time.min(),val_obs.time.max()))

    gr.timeseries(
                    val_obs[var],
                    cast=[clim_fc[var],hindcast[var]],
                    title='EC',
                    filename='ppEC',
                    clabs=['clim','EC']
                )
    gr.timeseries(
                    val_obs[var],
                    cast=[clim_fc[var],hindcast[var]],
                    lead_time=[28,35,42],
                    title='EC',
                    filename='ppEC',
                    clabs=['clim','EC']
                )

    crps_mod  = xs.crps_ensemble(val_obs[var],hindcast[var],dim=[])
    crps_clim = xs.crps_gaussian(
                                val_obs[var],
                                clim_mean[var],
                                clim_std[var],
                                dim=[]
                                )

    gr.skill_plot(crps_mod,crps_clim,title='EC',filename='ppEC')
    gr.qq_plot(val_obs[var],hindcast[var],y_axlabel='EC',filename='ppEC')
