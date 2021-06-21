import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import time as time_lib
import gridpp

from S2S.local_configuration import config
from S2S.data_handler        import BarentsWatch, ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.graphics.graphics as gr

import scripts.Henrik.create_domain_file

path_e = 't2m/'
long_name = 'absolute_t2m'
Archive().make_dir(config['VALID_DB']+path_e)

domainID = 'norwegian_coast'

#var      = 't2m'

#var1     = False
#var2     = False

var      = 'abs_wind'
var1     = 'u10'
var2     = 'v10'

t_start  = (2019,7,1)
t_end    = (2020,6,26)
# t_end    = (2020,2,3)


clim_t_start  = (1999,1,1)
clim_t_end    = (2021,1,4)

process_hindcast     = True
process_era          = True
make_time_series     = True

high_res             = False

# steps = pd.to_timedelta([7,14,23,30,37,44],'D')
steps = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')


    
    
print('\nProcess hindcast')
if process_hindcast:

    if var1 is False:
        
        print('\tLoad hindcast')
        hindcast = ECMWF_S2SH(high_res=high_res)\
                        .load(
                                var,
                                t_start,
                                t_end,
                                domainID,
                                download=False,
                                x_landmask=False
                            )[var]

        hindcast = hindcast.sortby(['time','step'])

        print('\tApply 7D running mean along lead time dimension')
        hindcast = hindcast.rolling(step=7,center=True)\
                                        .mean().dropna('step')

        print('\tKeep only lead times of interest')
        hindcast = hindcast.where(hindcast.step.isin(steps),drop=True)
    else:     
        print('\tLoad hindcast')
        hindcast_u = ECMWF_S2SH(high_res=high_res)\
                        .load(
                                var1,
                                t_start,
                                t_end,
                                domainID,
                                download=False,
                                x_landmask=False
                            )[var1]

        hindcast_u = hindcast_u.sortby(['time','step'])

        print('\tApply 7D running mean along lead time dimension')
        hindcast_u = hindcast_u.rolling(step=7,center=True)\
                                        .mean().dropna('step')

        print('\tKeep only lead times of interest')
        hindcast_u = hindcast_u.where(hindcast_u.step.isin(steps),drop=True)

        print('\tLoad hindcast')
        hindcast_v = ECMWF_S2SH(high_res=high_res)\
                        .load(
                                var2,
                                t_start,
                                t_end,
                                domainID,
                                download=False,
                                x_landmask=False
                            )[var2]

        hindcast_v = hindcast_v.sortby(['time','step'])

        print('\tApply 7D running mean along lead time dimension')
        hindcast_v = hindcast_v.rolling(step=7,center=True)\
                                        .mean().dropna('step')

        print('\tKeep only lead times of interest')
        hindcast_v = hindcast_v.where(hindcast_v.step.isin(steps),drop=True)

        print('\tCompute absolute wind')
        hindcast_u,hindcast_v = xr.align(hindcast_u,hindcast_v)

        hindcast = xr.apply_ufunc(
                        xh.absolute,hindcast_u,hindcast_v,
                        input_core_dims  = [[],[]],
                        output_core_dims = [[]],
                        vectorize=True,dask='parallelized'
                    )
        hindcast = hindcast.rename(var)

    try:
        hindcast = hindcast.drop('valid_time')
    except ValueError:
        pass

    try:
        hindcast = hindcast.drop('number')
    except ValueError:
        pass
    try:
        hindcast = hindcast.drop('surface')
    except ValueError:
        pass

    # store in absolute values
    hindcast.to_netcdf(config['VALID_DB']+path_e+long_name+'_hindcast.nc')

    print('\tCompute model climatology')
    mean,std = xh.c_climatology(hindcast)

    hindcast_a = (hindcast-mean)/std

    hindcast_a = hindcast_a.rename(var)

    # store hindcast anomlies
    hindcast_a.to_netcdf(config['VALID_DB']+path_e+long_name+\
                                '_anomalies_hindcast.nc')

hindcast = xr.open_dataset(config['VALID_DB']+path_e+long_name+\
                                '_hindcast.nc')[var]

hindcast_a = xr.open_dataset(config['VALID_DB']+path_e+long_name+\
                            '_anomalies_hindcast.nc')[var]

print('\nProcess ERA as observations')
if process_era:

    print('\tLoad ERA')
    era_u = ERA5(high_res=high_res)\
                .load(var1,clim_t_start,clim_t_end,domainID)[var1]

    era_v = ERA5(high_res=high_res)\
                .load(var2,clim_t_start,clim_t_end,domainID)[var2]

    era_u,era_v = xr.align(era_u,era_v)

    era = xr.apply_ufunc(
                    xh.absolute,era_u,era_v,
                    input_core_dims  = [[],[]],
                    output_core_dims = [[]],
                    vectorize=True,dask='parallelized'
                )

    era = era.rename(var)

    era = era.sortby(['time'])

    print('\tApply 7D running mean along time dimension')
    era = era.rolling(time=7,center=True).mean().dropna('time')

    print('\tOrganize ERA like hindcast')
    stacked_era = \
        xh.at_validation(era,hindcast.time+hindcast.step,ddays=1)

    print('\tCompute climatology')
    clim_mean,clim_std = xh.o_climatology(stacked_era)

    stacked_era_a = (stacked_era-clim_mean)/clim_std

    stacked_era   = stacked_era.drop('validation_time').rename(var)
    stacked_era_a = stacked_era_a.drop('validation_time').rename(var)
    clim_mean     = clim_mean.drop('validation_time').rename(var)
    clim_std      = clim_std.drop('validation_time').rename(var)

    stacked_era.to_netcdf(config['VALID_DB']+path_e+\
                                'absoulte_wind_era.nc')
    stacked_era_a.to_netcdf(config['VALID_DB']+path_e+\
                                'absolute_wind_anomalies_era.nc')
    clim_mean.to_netcdf(config['VALID_DB']+path_e+\
                                'absolute_wind_era_mean.nc')
    clim_std.to_netcdf(config['VALID_DB']+path_e+\
                                'absolute_wind_era_std.nc')

    print('\tGenerate random forecasts')
    random_fc_a = models.deterministic_gaussian_forecast(
                                                xr.full_like(clim_mean,0.),
                                                xr.full_like(clim_std,1.)
                                                )
    random_fc   = models.deterministic_gaussian_forecast(
                                                clim_mean,
                                                clim_std
                                                )

    random_fc_a.to_netcdf(config['VALID_DB']+path_e+\
                            'absolute_wind_random-forecast_anomalies.nc')[var]

    random_fc.to_netcdf(config['VALID_DB']+path_e+\
                                'absolute_wind_random-forecast.nc')[var]

if make_time_series:

    stacked_era = xr.open_dataset(config['VALID_DB']+path_e+\
                                            'absoulte_wind_era.nc')[var]

    stacked_era_a = xr.open_dataset(config['VALID_DB']+path_e+\
                                'absolute_wind_anomalies_era.nc')[var]

    clim_mean = xr.open_dataset(config['VALID_DB']+path_e+\
                                'absolute_wind_era_mean.nc')[var]

    clim_std = xr.open_dataset(config['VALID_DB']+path_e+\
                                'absolute_wind_era_std.nc')[var]

    random_fc = xr.open_dataset(config['VALID_DB']+path_e+\
                                'absolute_wind_random-forecast.nc')[var]

    random_fc_a = xr.open_dataset(config['VALID_DB']+path_e+\
                                'absolute_wind_random-forecast_anomalies.nc')[var]

    stacked_era = xh.assign_validation_time(
                    stacked_era.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    stacked_era_a = xh.assign_validation_time(
                    stacked_era_a.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    hindcast = xh.assign_validation_time(
                    hindcast.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    hindcast_a = xh.assign_validation_time(
                    hindcast_a.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    clim_mean = xh.assign_validation_time(
                    clim_mean.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    clim_std = xh.assign_validation_time(
                    clim_std.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    random_fc = xh.assign_validation_time(
                    random_fc.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    random_fc_a = xh.assign_validation_time(
                    random_fc_a.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    #[7, 14, 21, 28, 35, 42]
    gr.timeseries(
                    stacked_era,
                    cast=[random_fc],
                    title='EC',
                    filename='wind_random',
                    clabs=['random_fc'],
                    lead_time=[7, 14]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[random_fc,hindcast,hindcast_a*clim_std + clim_mean],
                    title='EC',
                    filename='wind_abs_all',
                    clabs=['random_fc','EC','EC_sb'],
                    lead_time=[7, 14]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[hindcast,hindcast_a*clim_std + clim_mean],
                    title='EC',
                    filename='wind_abs_fc',
                    clabs=['EC','EC_sb'],
                    lead_time=[7, 14]
                )

    gr.timeseries(
                    stacked_era_a,
                    cast=[random_fc_a,hindcast_a],
                    title='EC',
                    filename='wind_anom_all',
                    clabs=['random','EC_a'],
                    lead_time=[7, 14]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[random_fc],
                    title='EC',
                    filename='wind_random',
                    clabs=['random_fc'],
                    lead_time=[21, 28]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[random_fc,hindcast,hindcast_a*clim_std + clim_mean],
                    title='EC',
                    filename='wind_abs_all',
                    clabs=['random_fc','EC','EC_sb'],
                    lead_time=[21, 28]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[hindcast,hindcast_a*clim_std + clim_mean],
                    title='EC',
                    filename='wind_abs_fc',
                    clabs=['EC','EC_sb'],
                    lead_time=[21, 28]
                )

    gr.timeseries(
                    stacked_era_a,
                    cast=[random_fc_a,hindcast_a],
                    title='EC',
                    filename='wind_anom_all',
                    clabs=['random','EC_a'],
                    lead_time=[21, 28]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[random_fc],
                    title='EC',
                    filename='wind_random',
                    clabs=['random_fc'],
                    lead_time=[35,42]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[random_fc,hindcast,hindcast_a*clim_std + clim_mean],
                    title='EC',
                    filename='wind_abs_all',
                    clabs=['random_fc','EC','EC_sb'],
                    lead_time=[35,42]
                )

    gr.timeseries(
                    stacked_era,
                    cast=[hindcast,hindcast_a*clim_std + clim_mean],
                    title='EC',
                    filename='wind_abs_fc',
                    clabs=['EC','EC_sb'],
                    lead_time=[35,42]
                )

    gr.timeseries(
                    stacked_era_a,
                    cast=[random_fc_a,hindcast_a],
                    title='EC',
                    filename='wind_anom_all',
                    clabs=['random','EC_a'],
                    lead_time=[35,42]
                )
