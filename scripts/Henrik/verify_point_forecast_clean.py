import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs

from scripts.Henrik.data_handler import BarentsWatch, ERA5, ECMWF_S2SH, Archive
import scripts.Henrik.xarray_helpers as xh
from S2S.local_configuration import config
import scripts.Henrik.models as models

domainID = 'NVK'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2020,3,16)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

process_hindcast_and_training_data = False
train_models                       = False
process_point_observations         = False
verify                             = True

point_observations = BarentsWatch().load(['Hisdalen','Langøy S'])

if process_hindcast_and_training_data:
    ###################
    #### Load data ####
    ###################

    observations = ERA5().load(var,clim_t_start,clim_t_end,domainID)[var]-272.15

    hindcast     = ECMWF_S2SH().load(var,t_start,t_end,domainID)[var]-272.15

    ############################################################################
    #### Inter/extrapolate NaNs by nearest functioning gridpoint to the west ###
    #### Interpolate to locations of BW observations                         ###
    ############################################################################

    observations       = observations\
                            .interpolate_na(
                                            dim='lon',
                                            method='nearest',
                                            fill_value="extrapolate"
                                        )\
                                .interp(
                                        lon=point_observations.lon,
                                        lat=point_observations.lat,
                                        method='nearest'
                                    )

    hindcast           = hindcast\
                            .interpolate_na(
                                            dim='lon',
                                            method='nearest',
                                            fill_value="extrapolate"
                                        )\
                                .interp(
                                        lon=point_observations.lon,
                                        lat=point_observations.lat,
                                        method='nearest'
                                    )

    ##############################################
    #### Compute running means with 7D window ####
    ##############################################

    observations = observations\
                    .rolling(time=7,center=True).mean()\
                        .dropna('time')\

    hindcast     = hindcast\
                    .rolling(step=7,center=True).mean()\
                        .dropna('step')\

    #####################
    #### Climatology ####
    #####################
    obs_mean = observations.groupby('time.month').mean()
    obs_std  = observations.groupby('time.month').std()

    hc_mean  = hindcast.mean('member')
    hc_std   = hindcast.std('member')

    #####################
    #### Standardize ####
    #####################

    observations       = (
                            observations.groupby('time.month')-obs_mean
                        ).groupby('time.month')/obs_std

    hindcast           = (hindcast - hc_mean)/hc_std

    observations       = observations.to_dataset(name=var)
    hindcast           = hindcast.to_dataset(name=var)

    hindcast.to_netcdf(config['VALID_DB']+'/h_temp.nc')
    observations.to_netcdf(config['VALID_DB']+'/o_temp.nc')

hindcast = xr.open_dataset(config['VALID_DB']+'/h_temp.nc')
observations = xr.open_dataset(config['VALID_DB']+'/o_temp.nc')

# Match dimansions of hindcast and observations
hindcast,observations_val = xh.match_times(hindcast,observations)
# Observations at cast time stacked like hindcast
observations_init_time    = xh.stack_like(observations,hindcast)

if train_models:

    hindcast                   = hindcast.transpose(
                                            'member','time','step','location'
                                            )
    observations_val           = observations_val.transpose(
                                                        'time','step','location'
                                                        )

    p_mod = models.persistence(observations_init_time,observations_val,var)
    c_mod = models.combo(observations_init_time,hindcast,observations_val,var)

    p_mod.to_netcdf(config['VALID_DB']+'/p_mod_ERA_temp.nc')
    c_mod.to_netcdf(config['VALID_DB']+'/c_mod_ERA_temp.nc')

p_mod = xr.open_dataset(config['VALID_DB']+'/p_mod_ERA_temp.nc')
c_mod = xr.open_dataset(config['VALID_DB']+'/c_mod_ERA_temp.nc')

if verify:

    for loc in point_observations.location:

        point_observations = point_observations.sel(location=loc)
        point_observations = xh.keep_groups_of(
                                        point_observations,
                                        dim='time.month',
                                        members=11
                                    )
        clim_mean = point_observations.groupby('time.month').mean()
        clim_std  = point_observations.groupby('time.month').std()

        anomalies = (point_observations.groupby('time.month') - clim_mean)\
                            .groupby('time.month')/clim_std

        for step in hindcast.step:
            

exit()
#########################
#### Calculate skill ####
#########################

crps_fc   = xs.crps_ensemble(point_observations,point_hindcast,dim=[])
crps_fc_g = xs.crps_gaussian(point_observations,point_hindcast.mean('member'),point_hindcast.std('member'),dim=[])
crps_fc_g = xs.crps_gaussian(point_observations,point_hindcast.mean('member'),point_hindcast.std('member'),dim=[])
print(observations)
print(point_obs_std)
print(obs_std)
exit()
