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

high_res             = True

all_observations = BarentsWatch().load('all',no=200).sortby('time')
di = 5

for ii in range(di,len(all_observations.location),di):

    observations = all_observations.isel(location=slice(ii-di,ii))

    print('\nProcess hindcast')
    if process_hindcast:

        print('\tLoad hindcast')
        hindcast = ECMWF_S2SH(high_res=high_res)\
                        .load(var,t_start,t_end,domainID)[var]-272.15

        # print('\tExtrapolate land mask')
        # hindcast = xh.extrapolate_land_mask(hindcast)

        print('\tInterpolate hindcast to point locations')
        hindcast = xh.interp_to_loc(observations,hindcast)

        hindcast = hindcast.sortby(['time','step'])

        print('\tApply 7D running mean')
        hindcast = xh.running_mean(hindcast,window=7)

        print('\tAssign validation time')
        hindcast = xh.assign_validation_time(hindcast)

        xh.store_by_location(hindcast,'h_temp')

        process_observations = True

    hindcast = xh.load_by_location(observations.location,'h_temp')

    print('\nProcess observations')
    if process_observations:

        # standardize observations using full climatology
        observations_a = (
                            observations.groupby('time.month')\
                            -observations.groupby('time.month').mean(skipna=True)\
                        ).groupby('time.month')/\
                            observations.groupby('time.month').std(skipna=True)

        # calc climatology leaving current observation out
        clim_mean,clim_std = xh.o_climatology(observations[var],'time.month')

        val_obs   = xh.at_validation(observations,hindcast.validation_time)
        val_obs_a = xh.at_validation(observations_a,hindcast.validation_time)

        clim_mean = xh.at_validation(
                                clim_mean,
                                hindcast.validation_time,
                                ddays=3.5
                                )
        clim_std  = xh.at_validation(
                                clim_std,
                                hindcast.validation_time,
                                ddays=3.5
                                )

        xh.store_by_location(val_obs,'v_temp')
        xh.store_by_location(val_obs_a,'v_a_temp')

        xh.store_by_location(clim_mean,'clim_mean_temp')
        xh.store_by_location(clim_std,'clim_std_temp')

    val_obs   = xh.load_by_location(observations.location,'v_temp')
    val_obs_a = xh.load_by_location(observations.location,'v_a_temp')

    clim_mean = xh.load_by_location(observations.location,'clim_mean_temp')
    clim_std  = xh.load_by_location(observations.location,'clim_std_temp')

    print('\nCalibrate')
    if calibrate:

        mean,std = xh.c_climatology(hindcast,dim='validation_time.month')

        mean = xh.climatology_to_validation_time(mean,hindcast.time+hindcast.step)
        std  = xh.climatology_to_validation_time(std,hindcast.time+hindcast.step)

        hindcast = (hindcast-mean)/std

        xh.store_by_location(hindcast,'hp_temp')

        verify = True

    hindcast = xh.load_by_location(observations.location,'hp_temp')

    print('\nVerify')
    if verify:

        hindcast = hindcast.sel(time=slice(val_obs.time.min(),val_obs.time.max()))

        # # probabilistic skill
        # crps_mod  = xs.crps_ensemble(val_obs_a[var],hindcast[var],dim=[])
        # crps_clim = xs.crps_gaussian(
        #                             val_obs_a[var],
        #                             mu=0,
        #                             sig=1,
        #                             dim=[]
        #                             )
        #
        # gr.skill_plot(
        #                 crps_mod,
        #                 crps_clim,
        #                 title='EC',
        #                 filename='CRPSS_EC',
        #                 ylab='CRPSS'
        #             )

        # # deterministic skill
        # mae_mod  = xs.mae(val_obs_a[var],hindcast[var].mean('member'),dim=[])
        # mae_clim = xs.mae(
        #                   val_obs_a[var],
        #                   xr.full_like(hindcast[var].mean('member'),0),
        #                   dim=[]
        #                  )
        #
        # gr.skill_plot(
        #                 mae_mod,
        #                 mae_clim,
        #                 title='EC',
        #                 filename='MAESS_EC',
        #                 ylab='MAESS',
        #                 dim='validation_time.year'
        #             )

        # rmse_mod  = xs.rmse(val_obs_a[var],hindcast[var].mean('member'),dim=[])
        # rmse_clim = xs.rmse(
        #                   val_obs_a[var],
        #                   xr.full_like(hindcast[var].mean('member'),0),
        #                   dim=[]
        #                  )
        #
        # gr.skill_plot(
        #                 mae_mod,
        #                 mae_clim,
        #                 title='EC',
        #                 filename='RMSESS_EC',
        #                 ylab='RMSESS'
        #             )
        #
        # gr.qq_plot(val_obs_a[var],hindcast[var].mean('member'),
        #                             y_axlabel='EC',filename='EC')
        #
        # hindcast   = (hindcast*clim_std) + clim_mean
        # clim_fc    = models.clim_fc(clim_mean,clim_std)
        #
        #
        # gr.timeseries(
        #                 val_obs[var],
        #                 cast=[clim_fc[var],hindcast[var]],
        #                 title='EC',
        #                 filename='EC',
        #                 clabs=['clim','EC']
        #             )
        # gr.timeseries(
        #                 val_obs[var],
        #                 cast=[clim_fc[var],hindcast[var]],
        #                 lead_time=[28,35,42],
        #                 title='EC',
        #                 filename='EC',
        #                 clabs=['clim','EC']
        #             )
