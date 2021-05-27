import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import time as time_lib

from scripts.Henrik.data_handler import BarentsWatch, ERA5, ECMWF_S2SH, Archive
import scripts.Henrik.xarray_helpers as xh
from S2S.local_configuration import config
import scripts.Henrik.models as models
import scripts.Henrik.organize_barentzwatch as ob

domainID = 'NVK'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,18)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

process_hindcast_and_training_data   = False
train_models                         = False
do_modeling                          = False
plot                                 = False
skill                                = True

point_observations = BarentsWatch().load(['Hisdalen','Langøy S']).sortby('time')

print('Process model and training data')
if process_hindcast_and_training_data:
    ###################
    #### Load data ####
    ###################

    print('Process model and training data: load ERA5')
    observations = ERA5().load(var,clim_t_start,clim_t_end,domainID)[var]-272.15

    print('Process model and training data: load hindcast')
    hindcast     = ECMWF_S2SH().load(var,t_start,t_end,domainID)[var]-272.15

    ############################################################################
    #### Inter/extrapolate NaNs by nearest functioning gridpoint to the west ###
    #### Interpolate to locations of BW observations                         ###
    ############################################################################
    print('Process model and training data: interpolate observations to point locations')
    observations       = observations.sortby('time')\
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
    print('Process model and training data: interpolate hindcast to point locations')
    hindcast           = hindcast.sortby(['time','step'])\
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
    print('Process model and training data: compute 7D running means')
    observations = observations\
                    .rolling(time=7,center=True).mean()\
                        .dropna('time')

    hindcast     = hindcast\
                    .rolling(step=7,center=True).mean()\
                        .dropna('step')

    #####################
    #### Climatology ####
    #####################
    print('Process model and training data: calculate climatology')
    obs_clim_mean,obs_clim_std = xh.o_climatology(observations,'time.month')
    hin_clim_mean,hin_clim_std = xh.c_climatology(hindcast,'time.dayofyear')

    #####################
    #### Standardize ####
    #####################
    print('Process model and training data: compute anomalies')
    observations = (observations-obs_clim_mean)/obs_clim_std
    hindcast     = (hindcast-hin_clim_mean)/hin_clim_std

    observations = observations.to_dataset(name=var)
    hindcast     = hindcast.to_dataset(name=var)

    hindcast.to_netcdf(config['VALID_DB']+'/h_temp.nc')
    observations.to_netcdf(config['VALID_DB']+'/o_temp.nc')

hindcast = xr.open_dataset(config['VALID_DB']+'/h_temp.nc')
observations = xr.open_dataset(config['VALID_DB']+'/o_temp.nc')

print('Train models')
if train_models:

    observations.sortby('time')
    # Match dimensions of hindcast and observations
    hindcast,observations_val = xh.match_times(hindcast,observations)
    # Observations at cast time stacked like hindcast
    observations_init_time = observations\
                                .sel(time=slice(
                                            hindcast.time.min(),
                                            hindcast.time.max()
                                            )
                                ).broadcast_like(hindcast.mean('member'))

    print('Train models: fit models')
    p_mod = models.persistence(observations_init_time,observations_val,var)
    c_mod = models.combo(observations_init_time,hindcast,observations_val,var)

    p_mod.to_netcdf(config['VALID_DB']+'/p_mod_ERA_temp.nc')
    c_mod.to_netcdf(config['VALID_DB']+'/c_mod_ERA_temp.nc')

p_mod = xr.open_dataset(config['VALID_DB']+'/p_mod_ERA_temp.nc')
c_mod = xr.open_dataset(config['VALID_DB']+'/c_mod_ERA_temp.nc')

print('Model')
if do_modeling:

    for loc in point_observations.location:

        loc_str = str(loc.values)

        print('Model: pick location')
        p_obs = point_observations.sel(location=loc)

        print('Model: pick groups')
        p_obs = xh.keep_groups_of(
                                        p_obs,
                                        dim='time.month',
                                        members=11
                                    )

        print('Model: compute climatology')
        p_obs = p_obs.sortby('time')
        clim_mean,clim_std = xh.o_climatology(p_obs[var],'time.month')

        print('Model: calculate anomalies')
        p_obs_anom = (p_obs-clim_mean)/clim_std

        start = p_obs.time.min()
        end   = p_obs.time.max()

        print('Model: match times')
        hc                   = hindcast.sel(location=loc,time=slice(start,end))
        p_obs_init_time = p_obs_anom.reindex(
                                    indexers={'time':hc.time},
                                    method='pad',
                                    tolerance=pd.Timedelta(1,'W')
                                ).broadcast_like(hc.mean('member'))

        clim_mean = clim_mean.to_dataset(name=var)
        clim_std  = clim_std.to_dataset(name=var)

        clim_mean = clim_mean.resample(time='1D').nearest(tolerance='12H')
        clim_mean = xh.match_core(clim_mean,hc.time,hc.step)

        clim_std  = clim_std.resample(time='1D').nearest(tolerance='12H')
        clim_std  = xh.match_core(clim_std,hc.time,hc.step)

        p_obs = p_obs.resample(time='1D').nearest(tolerance='12H')
        p_obs = xh.match_core(p_obs,hc.time,hc.step)

        p_obs_anom = p_obs_anom.resample(time='1D').nearest(tolerance='12H')
        p_obs_anom = xh.match_core(p_obs_anom,hc.time,hc.step)

        start = p_obs_init_time.time.min()
        end   = p_obs_init_time.time.max()

        p_mod_l = p_mod.sel(location=loc,time=slice(start,end))
        c_mod_l = c_mod.sel(location=loc,time=slice(start,end))

        print('Model: model')
        model       = xr.merge(
                        [
                            (
                                clim_mean + hc[var]*clim_std
                            ).rename({var:'absolute'})\
                                .transpose('member','step','time'),
                            (
                                hc[var]
                            ).rename('anomalies')\
                                .transpose('member','step','time')
                        ], compat='override', join='outer'
                    )

        persistence = xr.merge(
                        [
                            (
                                clim_mean+\
                                (
                                p_mod_l.slope*p_obs_init_time[var]
                                )*clim_std
                            ).rename({var:'absolute'})\
                                .transpose('step','time'),
                            (
                                p_mod_l.slope*p_obs_init_time[var]
                            ).rename('anomalies')\
                                .transpose('step','time')
                        ], compat='override', join='outer'
                    )

        combo      = xr.merge(
                        [
                            (
                                clim_mean+\
                                (
                                c_mod_l.slope_obs*p_obs_init_time[var]+\
                                c_mod_l.slope_model*hc[var])*clim_std
                            ).rename({var:'absolute'})\
                                .transpose('member','step','time'),
                            (
                                c_mod_l.slope_obs+\
                                c_mod_l.slope_model*hc[var]
                            ).rename('anomalies')\
                                .transpose('member','step','time')
                        ], compat='override', join='outer'
                    )

        model.to_netcdf(config['VALID_DB']+'/computed_model_'+loc_str+'.nc')
        persistence.to_netcdf(config['VALID_DB']+'/computed_pers_'+loc_str+'.nc')
        combo.to_netcdf(config['VALID_DB']+'/computed_combo_'+loc_str+'.nc')

        xr.merge(
                        [
                            (
                                p_obs
                            ).rename({var:'absolute'})\
                                .transpose('step','time'),
                            (
                                p_obs_anom
                            ).rename({var:'anomalies'})\
                                .transpose('step','time')
                        ], compat='override', join='outer'
                    ).to_netcdf(config['VALID_DB']+\
                                '/processed_observations_'+loc_str+'.nc')

        xr.merge(
                        [
                            (
                                clim_mean
                            ).rename({var:'absolute'})\
                                .transpose('step','time'),
                            (
                                clim_std
                            ).rename({var:'anomalies'})\
                                .transpose('step','time')
                        ], compat='override', join='outer'
                    ).to_netcdf(config['VALID_DB']+\
                                '/processed_climatology_'+loc_str+'.nc')
if plot:

    for loc in point_observations.location:

        loc_str = str(loc.values)

        obs         = xr.open_dataset(config['VALID_DB']+\
                                    '/processed_observations_'+loc_str+'.nc')
        model       = xr.open_dataset(config['VALID_DB']+\
                                    '/computed_model_'+loc_str+'.nc')
        pers        = xr.open_dataset(config['VALID_DB']+\
                                    '/computed_pers_'+loc_str+'.nc')
        combo       = xr.open_dataset(config['VALID_DB']+\
                                    '/computed_combo_'+loc_str+'.nc')

        for lt in range(10,30):

            o   = obs.sel(step=pd.Timedelta(lt,'D')).dropna('time')
            m   = model.sel(step=pd.Timedelta(lt,'D')).dropna('time')
            p   = pers.sel(step=pd.Timedelta(lt,'D')).dropna('time')
            c   = combo.sel(step=pd.Timedelta(lt,'D')).dropna('time')

            plt.plot(o.time,o.absolute,'k',linewidth=3)
            plt.plot(m.time,m.mean('member').absolute,alpha=0.5,label='model')
            plt.plot(p.time,p.absolute,alpha=0.5,label='pers')
            plt.plot(c.time,c.mean('member').absolute,alpha=0.5,label='combo')
            plt.title(ob.name_from_loc(loc_str)+' Lead time: '+str(lt)+' days')
            plt.legend()
            plt.show()

if skill:

    for loc in point_observations.location:

        loc_str     = str(loc.values)

        obs         = xr.open_dataset(config['VALID_DB']+\
                                    '/processed_observations_'+loc_str+'.nc')
        model       = xr.open_dataset(config['VALID_DB']+\
                                    '/computed_model_'+loc_str+'.nc')
        pers        = xr.open_dataset(config['VALID_DB']+\
                                    '/computed_pers_'+loc_str+'.nc')
        combo       = xr.open_dataset(config['VALID_DB']+\
                                    '/computed_combo_'+loc_str+'.nc')
        clim        = xr.open_dataset(config['VALID_DB']+\
                                    '/processed_climatology_'+loc_str+'.nc')

        obs         = obs.isel(step=slice(10,-1))
        model       = model.isel(step=slice(10,-1))
        pers        = pers.isel(step=slice(10,-1))
        combo       = combo.isel(step=slice(10,-1))
        clim        = clim.isel(step=slice(10,-1))

        import scripts.Henrik.graphics as gr

        crps_model = xs.crps_ensemble(obs.anomalies,model.anomalies,dim=[])
        crps_combo = xs.crps_ensemble(obs.anomalies,combo.anomalies,dim=[])
        crps_clim  = xs.crps_gaussian(obs.anomalies,0,1,dim=[])

        gr.skill_plot(crps_model,crps_clim,title='EC',filename='ECcrpss')
        gr.skill_plot(crps_combo,crps_clim,title='COMBO',filename='COMBOcrpss')

        gr.qq_plot(obs.anomalies,model.anomalies,title='EC',filename='ECqq')
        gr.qq_plot(obs.anomalies,pers.anomalies, title='PERS',filename='PERSqq')
        gr.qq_plot(obs.anomalies,combo.anomalies,title='COMBO',filename='COMBOqq')
