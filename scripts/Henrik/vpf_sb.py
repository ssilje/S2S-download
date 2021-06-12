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

high_res             = True

hindcast = ECMWF_S2SH(high_res=high_res)\
                .load(
                        var,
                        t_start,
                        t_end,
                        domainID,
                        download=True,
                        x_landmask=True
                    )[var]-272.15

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

        # drop pre-assigned validation time coordinates
        hindcast = hindcast.drop('valid_time')

        print('\tInterpolate hindcast to point locations')
        hindcast = xh.interp_to_loc(observations,hindcast)

        hindcast = hindcast.sortby(['time','step'])

        print('\tApply 7D running mean along lead time dimension')
        hindcast = xh.running_mean(hindcast,window=7)

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
        stacked_obs = xh.at_validation(observations,hindcast.time+hindcast.step)

        print('\tCompute climatology')
        clim_mean,clim_std = xh.o_climatology(stacked_obs)

        stacked_obs_a = (stacked_obs-clim_mean)/clim_std

        xh.store_by_location(stacked_obs,'v_temp')
        xh.store_by_location(stacked_obs_a,'v_a_temp')

        xh.store_by_location(clim_mean,'clim_mean_temp')
        xh.store_by_location(clim_std,'clim_std_temp')

    # hindcast   = xh.load_by_location(observations.location,'h_temp')
    # hindcast_a = xh.load_by_location(observations.location,'ha_temp')
    #
    # stacked_obs   = xh.load_by_location(observations.location,'v_temp')
    # stacked_obs_a = xh.load_by_location(observations.location,'v_a_temp')
    #
    # clim_mean = xh.load_by_location(observations.location,'clim_mean_temp')
    # clim_std  = xh.load_by_location(observations.location,'clim_std_temp')
    #
    # clim_fc  = models.clim_fc(clim_mean,clim_std).to_array().rename(var)
    # hindcast = hindcast.sel(time=slice(stacked_obs.time.min(),stacked_obs.time.max()))
    #
    # hindcast      = xh.assign_validation_time(hindcast)
    # hindcast_a    = xh.assign_validation_time(hindcast_a)
    # clim_fc       = xh.assign_validation_time(clim_fc)
    # stacked_obs   = xh.assign_validation_time(stacked_obs)
    # stacked_obs_a = xh.assign_validation_time(stacked_obs_a)
    #
    # # gr.timeseries(
    # #                 stacked_obs[var],
    # #                 cast=[clim_fc,hindcast[var]],
    # #                 title='EC',
    # #                 filename='test',
    # #                 clabs=['clim','EC']
    # #             )
    #
    # gr.timeseries(
    #                 stacked_obs_a[var],
    #                 cast=[hindcast_a.to_array().rename(var)],
    #                 title='EC',
    #                 filename='test_anomaly',
    #                 clabs=['EC']
    #             )
    # exit()
