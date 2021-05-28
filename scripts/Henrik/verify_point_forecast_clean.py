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
import scripts.Henrik.graphics as gr
import scripts.Henrik.latex as latex

domainID = 'NVK'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2020,1,8)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

process_hindcast_and_training_data   = True
train_models                         = True
do_modeling                          = True
plot                                 = True
skill                                = True

point_observations = BarentsWatch().load(['Hisdalen','Langøy S']).sortby('time')

print('Process model and training data')
if process_hindcast_and_training_data:
    ###################
    #### Load data ####
    ###################

    print('Process model and training data: load ERA5')
    observations = ERA5(high_res=True).load(var,clim_t_start,clim_t_end,domainID)[var]-272.15

    print('Process model and training data: load hindcast')
    hindcast     = ECMWF_S2SH(high_res=True).load(var,t_start,t_end,domainID)[var]-272.15

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
    hin_clim_mean,hin_clim_std = xh.c_climatology(hindcast,'time.month')

    #####################
    #### Standardize ####
    #####################
    print('Process model and training data: compute anomalies')
    observations = (observations-obs_clim_mean)/obs_clim_std
    hindcast     = (hindcast-hin_clim_mean)/hin_clim_std

    # hindcast = xh.c_by_vt(hindcast)
    # observations.groupby('time.month').map(xh.standard)

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
    observations_init_time = observations.reindex(
                                            {'time':hindcast.time},
                                            method='pad',
                                            tolerance='1W'
                                ).broadcast_like(hindcast.mean('member'))

    print('Train models: fit models')
    p_mod = models.persistence(observations_init_time,observations_val,var,dim='time.month')
    c_mod = models.combo(observations_init_time,hindcast,observations_val,var,dim='time.month')

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

        clim_mean = clim_mean.resample(time='1D').nearest(tolerance='3.5D')
        clim_mean = xh.match_core(clim_mean,hc.time,hc.step)

        clim_std  = clim_std.resample(time='1D').nearest(tolerance='3.5D')
        clim_std  = xh.match_core(clim_std,hc.time,hc.step)

        p_obs = p_obs.resample(time='1D').nearest(tolerance='3.5D')
        p_obs = xh.match_core(p_obs,hc.time,hc.step)

        p_obs_anom = p_obs_anom.resample(time='1D').nearest(tolerance='3.5D')
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
                                p_mod_l.intercept +
                                p_mod_l.slope*p_obs_init_time[var]
                                )*clim_std
                            ).rename({var:'absolute'})\
                                .transpose('step','time'),
                            (
                                p_mod_l.intercept+p_mod_l.slope*p_obs_init_time[var]
                            ).rename('anomalies')\
                                .transpose('step','time')
                        ], compat='override', join='outer'
                    )

        combo      = xr.merge(
                        [
                            (
                                clim_mean+\
                                (
                                c_mod_l.intercept +
                                c_mod_l.slope_obs*p_obs_init_time[var]+\
                                c_mod_l.slope_model*hc[var])*clim_std
                            ).rename({var:'absolute'})\
                                .transpose('member','step','time'),
                            (
                                c_mod_l.intercept+
                                c_mod_l.slope_obs*p_obs_init_time[var]+\
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
                            ).rename({var:'mean_value'})\
                                .transpose('step','time'),
                            (
                                clim_std
                            ).rename({var:'std_value'})\
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

        # for time in obs.time:
        #     o = obs.sel(time=time).sortby('step')
        #     m = model.sel(time=time).sortby('step')
        #     c = combo.sel(time=time).sortby('step')
        #     plt.plot(time+o.step,o.absolute,'o',label='observasjoner')
        #     plt.plot(time+m.step,m.absolute.mean('member'),label='melding')
        #     plt.plot(time+c.step,c.absolute.mean('member'),label='cmelding')
        #     plt.legend()
        #     plt.show()
        # plt.close()
        # exit()

        latex.set_style(style='white')
        fig,axes = plt.subplots(3,1,
                        figsize=latex.set_size(width='thesis',
                            subplots=(1,1))
                        )

        for n,lt in enumerate([25,32,40]):
            ax = axes[n]
            o   = obs.sel(step=pd.Timedelta(lt,'D')).groupby('time.dayofyear').mean(skipna=True)
            m   = model.sel(step=pd.Timedelta(lt,'D')).groupby('time.dayofyear').mean(skipna=True)
            p   = pers.sel(step=pd.Timedelta(lt,'D')).groupby('time.dayofyear').mean(skipna=True)
            c   = combo.sel(step=pd.Timedelta(lt,'D')).groupby('time.dayofyear').mean(skipna=True)


            # plt.plot(m.time,m.mean('member').absolute,alpha=0.5,label='model')
            # plt.plot(p.time,p.absolute,alpha=0.5,label='pers')
            ax.plot(c.dayofyear,c.mean('member').absolute,alpha=0.7,
                            label='combo',linewidth=0.5,zorder=30)
            ax.fill_between(c.dayofyear,
                                c.max('member').absolute,
                                c.min('member').absolute,alpha=0.3,zorder=30)

            ax.plot(m.dayofyear,m.mean('member').absolute,alpha=0.7,label='ec',linewidth=0.5)
            ax.fill_between(m.dayofyear,
                                m.max('member').absolute,
                                m.min('member').absolute,alpha=0.3)

            ax.plot(o.dayofyear,o.absolute,'k',linewidth=0.5,label='barentswatch')
            # plt.plot(c.time,c.min('member').absolute,'red',alpha=0.5)
            # plt.plot(c.time,c.max('member').absolute,'red',alpha=0.5)
            ax.set_title(ob.name_from_loc(loc_str)+' Lead time: '+str(lt)+' days')
            ax.legend()
            ax.set_xlabel('day of year')
            ax.set_xlabel('SST')
        gr.save_fig(fig,ob.name_from_loc(loc_str)+'_timeseries_mean')

        latex.set_style(style='white')
        fig,axes = plt.subplots(3,1,
                        figsize=latex.set_size(width='thesis',
                            subplots=(1,1))
                        )

        for n,lt in enumerate([25,32,40]):
            ax = axes[n]
            o   = obs.sel(step=pd.Timedelta(lt,'D'))
            m   = model.sel(step=pd.Timedelta(lt,'D'))
            p   = pers.sel(step=pd.Timedelta(lt,'D'))
            c   = combo.sel(step=pd.Timedelta(lt,'D'))


            # plt.plot(m.time,m.mean('member').absolute,alpha=0.5,label='model')
            # plt.plot(p.time,p.absolute,alpha=0.5,label='pers')
            ax.plot(c.time,c.mean('member').absolute,alpha=0.7,
                            label='combo',linewidth=0.5,zorder=30)
            ax.fill_between(c.time,
                                c.max('member').absolute,
                                c.min('member').absolute,alpha=0.3,zorder=30)

            ax.plot(m.time,m.mean('member').absolute,alpha=0.7,label='ec',linewidth=0.5)
            ax.fill_between(m.time,
                                m.max('member').absolute,
                                m.min('member').absolute,alpha=0.3)

            ax.plot(o.time,o.absolute,'k',linewidth=0.5,label='barentswatch')
            # plt.plot(c.time,c.min('member').absolute,'red',alpha=0.5)
            # plt.plot(c.time,c.max('member').absolute,'red',alpha=0.5)
            ax.set_title(ob.name_from_loc(loc_str)+' Lead time: '+str(lt)+' days')
            ax.legend()
            ax.set_xlabel('time')
            ax.set_xlabel('SST')
        gr.save_fig(fig,ob.name_from_loc(loc_str)+'_timeseries')

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

        crps_model = xs.crps_ensemble(obs.absolute,model.absolute,dim=[])
        crps_combo = xs.crps_ensemble(obs.absolute,combo.absolute,dim=[])
        crps_clim  = xs.crps_gaussian(obs.absolute,
                                        clim.mean_value,
                                        clim.std_value,
                                        dim=[]
                                    )
        crps_mix   = xs.crps_gaussian(obs.absolute,
                                        combo.absolute.mean('member'),
                                        clim.std_value,
                                        dim=[]
                                    )

        gr.skill_plot(crps_model,crps_clim,title='EC',filename='ECcrpss')
        gr.skill_plot(crps_combo,crps_clim,title='COMBO',filename='COMBOcrpss')
        gr.skill_plot(crps_mix,crps_clim,title='COMBO-MIX',filename='COMBO_MIXcrpss')

        po = point_observations.sel(location=loc)[var]
        clim_mean,clim_std = xh.o_climatology(po,'time.month')
        po = (po-clim_mean)/clim_std
        o  = observations.sel(location=loc)[var]
        o  = o.reindex({'time':po.time},method='nearest',tolerance='1D')
        gr.qq_plot(o,po,dim='time.month',
                                x_axlabel='ERA5',
                                y_axlabel='Barentswatch',
                                filename='qq_plot_observations'
                                )
        gr.qq_plot(obs.anomalies,model.anomalies,title='EC',filename='ECqq')
        gr.qq_plot(obs.anomalies,pers.anomalies, title='PERS',filename='PERSqq')
        gr.qq_plot(obs.anomalies,combo.anomalies,title='COMBO',filename='COMBOqq')

        gr.qq_plot(obs.anomalies,model.anomalies.mean('member'),title='EC',filename='ECqq_mean')
        # gr.qq_plot(obs.anomalies,pers.anomalies, title='PERS',filename='PERSqq_mean')
        gr.qq_plot(obs.anomalies,combo.anomalies.mean('member'),title='COMBO',filename='COMBOqq_mean')
