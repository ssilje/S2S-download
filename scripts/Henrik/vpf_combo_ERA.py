
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

domainID = 'NVK'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

process_hindcast     = False
process_observations = False
get_training_data    = False
training_data        = 'ERA'
calibrate            = False
fit_model            = True
assemble_model       = False
verify               = True

high_res             = False

observations = BarentsWatch().load(['Hisdalen','Langøy S']).sortby('time')

print('Process hindcast')
if process_hindcast:

    print('\tLoad hindcast')
    hindcast = ECMWF_S2SH(high_res=high_res)\
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
                                method='linear'
                            )
    print('\tApply 7D running mean')
    hindcast = hindcast\
                    .rolling(step=7,center=True).mean()\
                        .dropna('step')

    print('\tAssign validation time')
    hindcast = xh.assign_validation_time(hindcast)

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

print('Process training data')
if training_data=='ERA':
    training_filename = '/ERA_tdata.nc'

if training_data=='BW':
    training_filename = '/BW_tdata.nc'

if get_training_data:

    if training_data=='ERA':

        tolerance = 1

        print('\tLoad ERA')
        tdata = ERA5(high_res=high_res)\
                    .load(var,clim_t_start,clim_t_end,domainID)[var]-272.15

        print('\tInterpolate ERA to point locations')
        tdata = tdata.sortby('time')\
                                .interpolate_na(
                                                dim='lon',
                                                method='nearest',
                                                fill_value="extrapolate"
                                            )\
                                    .interp(
                                            lon=observations.lon,
                                            lat=observations.lat,
                                            method='linear'
                                        )

        print('\tCompute 7D running means on ERA')
        tdata = tdata.rolling(time=7,center=True).mean().dropna('time')

    if training_data=='BW':

        tolerance = 2

        print('\tGet BW')
        tdata = observations[var]


    print('\tCompute anomlies')
    tdata = (tdata.groupby('time.month')\
                -tdata.groupby('time.month').mean(skipna=True))\
                    .groupby('time.month')\
                        /tdata.groupby('time.month').std(skipna=True)

    print('\tMatch to hindcast')
    tdatav = xh.at_validation(tdata,hindcast.validation_time,ddays=tolerance)
    tdatai = tdata.reindex(
                        {'time':hindcast.time},method='pad',tolerance='1W'
                    ).broadcast_like(tdatav)

    tdata = xr.merge([tdatav.rename('v'),tdatai.rename('i')],compat='override')

    tdata.drop('month').to_netcdf(config['VALID_DB']+training_filename)

    fit_model = True

tdata = xr.open_dataset(config['VALID_DB']+training_filename)

print('Fit model')
model_filename = '/c_mod_'+training_data+'_temp.nc'

if fit_model:

    print('\tFitting COMBO model')
    c_mod = models.combo(
                        tdata.i,
                        hindcast[var],
                        tdata.v,
                        dim='validation_time.month'
                    )

    c_mod = xh.assign_validation_time(c_mod.drop('validation_time'))

    c_mod.to_netcdf(config['VALID_DB']+model_filename)

    assemble_model = True

c_mod = xr.open_dataset(config['VALID_DB']+model_filename)

print('Assemble model')
combo_filename = '/combo_'+training_data+'.nc'

if assemble_model:

    print('\tGet observation at initialization time')
    obs_t0 = observations[var].reindex(
                        {'time':hindcast.time},method='pad',tolerance='1W'
                    ).broadcast_like(hindcast[var])

    print('\tModel')
    combo  =    c_mod.intercept\
              + c_mod.model_t*hindcast[var]\
              + c_mod.obs_t0*obs_t0

    print('\tRe-adjust model spread')
    mean,std = xh.c_climatology(combo,'validation_time.month')
    mean     = xh.climatology_to_validation_time(mean,combo.validation_time)
    std      = xh.climatology_to_validation_time(std,combo.validation_time)

    mean = mean.rename(var)
    std  = std.rename(var)

    # combo = (combo-mean)/std + mean

    combo = combo.rename(var)

    combo.to_netcdf(config['VALID_DB']+combo_filename)

    verify = True

combo = xr.open_dataset(config['VALID_DB']+combo_filename)

print('Verify')
if verify:

    hindcast   = hindcast*clim_std + clim_mean
    combo      = combo*clim_std + clim_mean
    clim_fc    = models.clim_fc(clim_mean,clim_std)
    hindcast   = hindcast.sel(time=slice(val_obs.time.min(),val_obs.time.max()))

    gr.timeseries(
                    val_obs[var],
                    cast=[clim_fc[var],hindcast[var],combo[var]],
                    title='EC',
                    filename='comboEC',
                    clabs=['clim','EC','COMBO']
                )
    gr.timeseries(
                    val_obs[var],
                    cast=[clim_fc[var],hindcast[var],combo[var]],
                    lead_time=[28,35,42],
                    title='EC',
                    filename='comboEC',
                    clabs=['clim','EC','COMBO']
                )

    crps_mod  = xs.crps_ensemble(val_obs[var],combo[var],dim=[])
    crps_clim = xs.crps_gaussian(
                                val_obs[var],
                                clim_mean[var],
                                clim_std[var],
                                dim=[]
                                )

    gr.skill_plot(crps_mod,crps_clim,title='COMBO',filename='comboEC')
    gr.qq_plot(val_obs[var],hindcast[var],y_axlabel='COMBO',filename='comboEC')
