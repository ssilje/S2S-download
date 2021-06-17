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
# t_end    = (2020,2,3)


clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

process_hindcast     = True
process_observations = False
process_era          = False
persistence_model    = True

high_res             = False

# steps = pd.to_timedelta([4,9,11,16,18,23,25,30,32,37,39,44],'D')
steps = pd.to_timedelta([9,16,23,30,37,44],'D')

all_observations = BarentsWatch().load('all',no=350).sortby('time')
di = 5

for ii in range(di,len(all_observations.location),di):

    xh.print_progress(ii,len(all_observations.location))

    observations = all_observations.isel(location=slice(ii-di,ii))

    print('\nProcess hindcast')
    if process_hindcast:

        print('\tLoad hindcast')
        hindcast = ECMWF_S2SH(high_res=high_res)\
                        .load(
                                var,
                                t_start,
                                t_end,
                                domainID,
                                download=False,
                                x_landmask=True
                            )[var]-272.15
        gr.point_map(all_observations,poi='Hisdalen')
        exit()
        # drop pre-assigned validation time coordinates
        hindcast = hindcast.drop('valid_time')

        print('\tInterpolate hindcast to point locations')
        hindcast = xh.interp_to_loc(observations,hindcast)

        hindcast = hindcast.sortby(['time','step'])

        print('\tApply 7D running mean along lead time dimension')
        hindcast = xh.running_mean(hindcast,window=7)

        hindcast = hindcast.where(hindcast.step.isin(steps),drop=True)

        try:
            hindcast = hindcast.drop('number')
        except ValueError:
            pass
        try:
            hindcast = hindcast.drop('surface')
        except ValueError:
            pass

        # store hindcast in absolute values
        xh.store_by_location(hindcast,'h_temp')

        print('\tCompute model climatology')
        mean,std = xh.c_climatology(hindcast)

        hindcast_a = (hindcast-mean)/std

        # store hindcast anomlies
        xh.store_by_location(hindcast_a.rename(var),'ha_temp')

    hindcast = xh.load_by_location(observations.location,'h_temp')

    print('\nProcess observations')
    if process_observations:

        print('\tOrganize observations like hindcast')
        stacked_obs = \
            xh.at_validation(observations,hindcast.time+hindcast.step,ddays=1)

        print('\tCompute climatology')
        clim_mean,clim_std = xh.o_climatology(stacked_obs)

        stacked_obs_a = (stacked_obs-clim_mean)/clim_std

        xh.store_by_location(stacked_obs.drop('validation_time'),'v_temp')
        xh.store_by_location(stacked_obs_a.drop('validation_time'),'v_a_temp')

        xh.store_by_location(clim_mean.drop('validation_time'),'clim_mean_temp')
        xh.store_by_location(clim_std.drop('validation_time'),'clim_std_temp')

    print('\nProcess ERA as observations')
    if process_era:

        print('\tLoad ERA')
        era = ERA5(high_res=high_res)\
                    .load(var,clim_t_start,clim_t_end,domainID)[var]-272.15

        print('\tInterpolate ERA to point locations')
        era = xh.interp_to_loc(observations,era)

        era = era.sortby(['time'])

        print('\tApply 7D running mean along time dimension')
        era = xh.running_mean_o(era,window=7)

        xh.store_by_location(era,'era_temp')

        print('\tOrganize ERA like hindcast')
        stacked_obs = \
            xh.at_validation(era,hindcast.time+hindcast.step,ddays=1)

        print('\tCompute climatology')
        clim_mean,clim_std = xh.o_climatology(stacked_obs)

        stacked_obs_a = (stacked_obs-clim_mean)/clim_std

        stacked_obs   = stacked_obs.drop('validation_time').rename(var)
        stacked_obs_a = stacked_obs_a.drop('validation_time').rename(var)
        clim_mean     = clim_mean.drop('validation_time').rename(var)
        clim_std      = clim_std.drop('validation_time').rename(var)

        xh.store_by_location(stacked_obs,'era_v_temp')
        xh.store_by_location(stacked_obs_a,'era_v_a_temp')

        xh.store_by_location(clim_mean,'era_clim_mean_temp')
        xh.store_by_location(clim_std,'era_clim_std_temp')

    print('\tPrepare persistence model initialization data')
    if persistence_model:

        # BW
        o = observations
        cm,cs = xh.o_climatology(o)
        o_a   = (o-cm)/cs
        o_a = o_a.reindex(
                {'time':hindcast.time},
                method='pad',
                tolerance='1W'
                ).broadcast_like(hindcast.mean('member'))
        xh.store_by_location(o_a,'obs_init_temp')

        # ERA
        o     = xh.load_by_location(observations.location,'era_temp')[var]\
                                                                .sortby('time')
        cm,cs = xh.o_climatology(o)

        cm = cm.groupby('time').mean(skipna=True)
        cs = cs.groupby('time').mean(skipna=True)

        o_a   = (o-cm)/cs
        o_a = o_a.reindex(
                {'time':hindcast.time},
                method='pad',
                tolerance='D'
                ).broadcast_like(hindcast.mean('member'))
        xh.store_by_location(o_a,'era_init_temp')
