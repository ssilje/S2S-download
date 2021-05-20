import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs

import scripts.Henrik.handle_domain as hd
from S2S.local_configuration import config
import scripts.Henrik.handle_datetime as dt
import scripts.Henrik.figure_receipts as fr
from .acc import ACC
from .create_domain_file import make_dir
from scripts.Henrik.data_handler import BarentsWatch, ERA5, ECMWF_S2SH
from scripts.Henrik.prepare_data import get_observations, get_hindcast
import scripts.Henrik.build_a_bear as build_a_bear
import scripts.Henrik.compute_reliability as cr
import scripts.Henrik.organize_barentzwatch as org_bw

################################################################################

def stack_times(ds1,ds2,dim_name='step'):
    """
    Stack ds1 along dim to match ds2.

    args:
        ds1:        xarray.Dataset with time dimension
        ds2:        xarray.Dataset with time and timedelta dimension
        dim_name:   must be step

    returns:
        stacked ds1: xarray.Dataset
    """

    out = []
    # print(ds2.time + )
    for t in ds2.time:

        validation_time = t + ds2.sel(time=t).step
        out.append(ds1.sel(time=validation_time))

    return xr.concat(out,dim_name)

def cap_time_to_match(cast,obs):
    """
    Cap of time to match the shortest ends of cast and obs.

    args:
        cast: xarray.Dataset with time and step dimension
        obs: xarray.Dataset with time dimension

    returns
        cast: xarray.Dataset
        obs: xarray.Dataset
    """

    if cast.time[0]-obs.time[0]<pd.Timdelta(0): # obs starts later
        cast = cast.sel(time=slice(find_matching(cast.time,obs.time[0],[1]),-1))
    else:
        obs = obs.sel(time=slice(find_matching(obs.time,cast.time[0],[1]),-1))

    if cast.time[-1]+cast.step[-1]-obs.time[-1]<pd.Timdelta(0): # cast ends first
        obs = obs.sel(time=slice(0,find_matching(obs.time,
                                                    cast.time[-1]+cast.step[-1],
                                                    [-1])))
    else:
        cast = cast.sel(time=slice(0,find_matching(cast.time,
                                                    obs.time[-1]-cast.step[-1],
                                                    [-1])))
    return cast,obs

def find_matching(array,*args):
    """
    Returns indices of closest (matching) values in array to *args
    """
    out = []

    if hasattr(args[-1],'__iter__'):
        targets = args[:-1]
        side    = np.array(args[-1])
    else:
        targets = args
        side    = np.zeros_like(len(args))

    for n,target in enumerate(targets):

        arr = np.array(array - target)

        if side[n] < 0:
            arr = arr[arr<pd.Timedelta(0)]
        elif side[n] > 0:
            arr = arr[arr>pd.Timedelta(0)]
        else:
            pass

        out.append(np.argmin(np.abs(arr)))

    return tuple(out)

def cap_cast_time(cast,obs):
    """
    Cap of time of cast to not exceed obs.

    args:
        cast: xarray.Dataset with time and step dimension
        obs: xarray.Dataset with time dimension

    returns
        cast: xarray.Dataset
    """
    cast = cast.sortby('time')
    start,end = find_matching(
                            cast.time,
                            obs.time[0],
                            obs.time[-1]-cast.step[-1],
                            [1,-1]
                        )
    return cast.isel(time=slice(start,end))
################################################################################

def match_times(cast,obs):
    """
    Match cast to obs and stack obs like cast

    args:
        cast: xarray.Dataset with time and step dimension
        obs: xarray.Dataset with time dimension

    returns
        obs: xarray.Dataset with time and step dimension
    """
    step  = cast.step.sortby('step').to_pandas()
    obs   = obs.sortby('time')

    start = obs.time[0].to_pandas()
    end   = obs.time[-1].to_pandas()-step.max()

    cast = cast.sortby('time').sel(
                                time=slice(
                                    start,
                                    end
                                )
                            )
    out_obs  = []
    obs = obs.reindex(
            time=pd.date_range(
                start=start,
                end=end,
                freq=step[2]-step[1]
            )
        )

    for time in cast.time:
        out_obs.append(obs.sel(time=time+cast.step).drop('time'))
    return cast,xr.concat(out_obs,cast.time)

def stack_climatology(clim,cast):
    """
    Stack climatology like cast

    args:
        cast: xarray.Dataset with time and step dimension
        clim: xarray.Dataset with day/month-ofyear dimension

    returns
        clim: xarray.Dataset with time and step dimension
    """
    clim   = clim.sortby('dayofyear')
    clim   = clim.reindex(dayofyear=np.arange(1,367,1,dtype='int'))
    t_out  = []
    for time in cast.time:
        s_out = []
        for step in cast.sel(time=time).step:
            date = pd.to_datetime(
                        xr.DataArray(
                            time.variable+step.variable
                        ).to_pandas()
                    )

            s_out.append(clim.sel(dayofyear=date.day_of_year))
        t_out.append(xr.concat(s_out,cast.step))
    return xr.concat(t_out,cast.time)


domainID = 'NVK'

t_start  = (2020,1,23)
t_end    = (2021,1,14)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

location = 'Hisdalen'

from .data_handler import BarentsWatch

if 1:
    ###################
    #### Load data ####
    ###################

    # org_bw.organize_files()
    point_observations = BarentsWatch().load(['Hisdalen','Langøy S'])

    #############################################################################
    #### Inter/extrapolate NaNs by nearest functioning gridpoint to the west ####
    #### Interpolate to locations of BW observations                         ####
    #############################################################################
    observations       = get_observations(domainID,clim_t_start,clim_t_end)\
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

    hindcast           = get_hindcast(domainID,t_start,t_end)\
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

    ############################################
    #### Compute running means of 7D window ####
    ############################################
    observations = observations.sst\
                    .rolling(time=7,center=True).mean()\
                        .dropna('time')\
                            .to_dataset(name='sst')


    hindcast     = hindcast.sst\
                    .rolling(step=7,center=True).mean()\
                        .dropna('step')\
                            .to_dataset(name='sst')
    #####################
    #### Climatology ####
    #####################
    point_obs_mean = point_observations.groupby('time.dayofyear').mean(skipna=True) # not necissary to store this variable
    point_obs_std  = point_observations.groupby('time.dayofyear').std(skipna=True)

    obs_mean = observations.groupby('time.dayofyear').mean() # not necissary to store this variable
    obs_std  = observations.groupby('time.dayofyear').std()

    hc_mean  = hindcast.groupby('time.dayofyear').mean().mean('member') # not necissary to store this variable

    obs_std  = obs_std.sst\
                .pad(dayofyear=7,mode='wrap')\
                    .rolling(dayofyear=7,center=True).mean()\
                        .isel(dayofyear=slice(7,-7))\
                            .to_dataset(name='sst')

    point_obs_std  = point_obs_std.sst\
                        .pad(dayofyear=7,mode='wrap')\
                            .rolling(dayofyear=7,center=True).mean()\
                                .isel(dayofyear=slice(7,-7))\
                                    .to_dataset(name='sst')

    ###########################
    #### Compute anomalies ####
    ###########################
    point_observations = point_observations.groupby('time.dayofyear') - point_obs_mean
    observations       = observations.groupby('time.dayofyear') - obs_mean
    hindcast           = hindcast.groupby('time.dayofyear') - hc_mean

    print(point_observations,observations,hindcast)
    ##############################################################
    #### Match times and stack observations to match hindcast ####
    ##############################################################
    point_hindcast,point_observations = match_times(hindcast,point_observations)
    hindcast,observations             = match_times(hindcast,observations)

    point_hindcast     = point_hindcast.transpose('member','time','step','loc')
    point_observations = point_observations.transpose('time','step','loc')

    hindcast           = hindcast.transpose('member','time','step','loc')
    observations       = observations.transpose('time','step','loc')

    point_obs_std = stack_climatology(point_obs_std,point_hindcast)
    obs_std       = stack_climatology(obs_std,hindcast)

    ##########################
    #### Temporarly store ####
    ##########################
    point_hindcast.to_netcdf(config['VALID_DB']+'/ph_temp.nc')
    point_observations.to_netcdf(config['VALID_DB']+'/oh_temp.nc')

    hindcast.to_netcdf(config['VALID_DB']+'/h_temp.nc')
    observations.to_netcdf(config['VALID_DB']+'/o_temp.nc')

    obs_std.to_netcdf(config['VALID_DB']+'/ostd_temp.nc')
    point_obs_std.to_netcdf(config['VALID_DB']+'/postd_temp.nc')

point_hindcast = xr.open_dataset(config['VALID_DB']+'/ph_temp.nc')
point_observations = xr.open_dataset(config['VALID_DB']+'/oh_temp.nc')

hindcast = xr.open_dataset(config['VALID_DB']+'/h_temp.nc')
observations = xr.open_dataset(config['VALID_DB']+'/o_temp.nc')

obs_std = xr.open_dataset(config['VALID_DB']+'/ostd_temp.nc')
point_obs_std = xr.open_dataset(config['VALID_DB']+'/postd_temp.nc')

#########################
#### Calculate skill ####
#########################

crps_fc   = xs.crps_ensemble(point_observations,point_hindcast,dim=[])
crps_fc_g = xs.crps_gaussian(point_observations,point_hindcast.mean('member'),point_hindcast.std('member'),dim=[])
print(observations)
print(point_obs_std)
print(obs_std)
exit()
################################################################################
################################################################################
################################################################################


validation_time       = hindcast.step + hindcast.time
point_validation_time = point_hindcast.step + point_hindcast.time

obs_anom_stck = obs_anom_stck.expand_dims('validation_time')\
                        .assign_coords(validation_time = validation_time)

hc_anom_stck = hc_anom.expand_dims('validation_time')\
                        .assign_coords(validation_time = validation_time)

#####################################
#### Compute anomalies and modes ####
#####################################
obs_mean = observations.groupby('time.dayofyear').mean() # not necissary to store this variable
obs_std  = observations.groupby('time.dayofyear').std()
obs_anom = observations.groupby('time.dayofyear') - obs_mean

hc_mean  = hindcast.groupby('time.dayofyear').mean().mean('ensemble_member') # not necissary to store this variable
hc_anom  = hindcast.groupby('time.dayofyear') - hc_mean

#### Climatological running mean ####
obs_std  = obs_std.sst\
            .pad(dayofyear=7,mode='wrap')\
                .rolling(dayofyear=7,center=True).mean()\
                    .isel(dayofyear=slice(7,-7))\
                        .to_dataset(name='sst')

#####################################
#### Get dimensions of hindcasts ####
#####################################
lon       = np.array(hc_anom.lon)
lat       = np.array(hc_anom.lat)
step0     = pd.to_timedelta(hc_anom.step)[0].days
time      = np.array([ct+np.array(hc_anom.step)\
                                    for ct in np.array(hc_anom.time)])

#################
#### Weights ####
#################
weights   = np.cos(
                np.deg2rad(
                    np.array([
                        np.array([
                            np.array([
                                np.meshgrid(hc_anom.lat,hc_anom.lon)[0]
                            ]*hc_anom.sizes['time'])
                        ]*hc_anom.sizes['step'])
                    ]*hc_anom.sizes['ensemble_member'])
                )
            )

############################################################
#### Stack to match hindcast and convert to numpy.array ####
############################################################
sst_obs  = np.stack([obs_anom.sst.sel(time=t) for t in time],axis=1)
std_obs  = np.stack([obs_std.sst.sel(dayofyear=pd.to_datetime(t).dayofyear)\
                                    for t in time],axis=1)
sst_hc   = np.array(hc_anom.sst)

###########################################################
#### Stack to match, but do not convert to numpy array #### Consider this for remaining operations using xarray and xskillscore
###########################################################
out = []

for step in hc_anom.step:

    validation_time = hc_anom.sel(ensemble_member=1).sel(step=step).time
    out.append(obs_anom.sel(time=validation_time).drop('ensemble_member'))

obs_anom_stck = xr.concat(out,'step').transpose('step','time','lon','lat')

##############################################
#### Include validation time as dimension ####
##############################################
validation_time=(obs_anom_stck.step + obs_anom_stck.time)

obs_anom_stck = obs_anom_stck.expand_dims('validation_time')\
                        .assign_coords(validation_time = validation_time)

hc_anom_stck = hc_anom.expand_dims('validation_time')\
                        .assign_coords(validation_time = validation_time)

#####################
#### Reliability ####
#####################
thresholds =\
    obs_anom.sst.mean('lon').mean('lat').groupby('time.season')\
            .quantile([0.75,0.25])

rel,sns,stp = cr.grouped_reliability_of_forecast(
                    forecast=hc_anom_stck.sst.mean('lon').mean('lat'),
                    observations=obs_anom_stck.sst.mean('lon').mean('lat'),
                    thresholds=thresholds,
                    direction=1,
                    dim='time',
                    group='validation_time.season'
                    )

fr.reliability_plot(rel,sns,pd.to_timedelta(stp))
exit()
#############################
#### Compute CRPS and SS ####
#############################
crps_fc   = ps.crps_ensemble(
                observations=sst_obs,
                forecasts=np.moveaxis(sst_hc,0,-1)
                )

crps_clim      = ps.crps_gaussian(x=sst_obs, mu=0, sig=std_obs) # does not weight with cos(lat), might be problematic in next operation crps_SS
crps_clim_stat = ps.crps_ensemble(
                        observations=sst_obs,
                        forecasts=np.zeros_like(sst_obs)
                        )

crps_ss   = 1 - np.nanmean(crps_fc,axis=(-1,-2))/np.nanmean(crps_clim,axis=(-1,-2))
crps_ss_stat   = 1 - np.nanmean(crps_fc,axis=(-1,-2))/np.nanmean(crps_clim_stat,axis=(-1,-2))

#############
#### ACC ####
#############
acc = ACC(sst_hc,sst_obs,weights)

#################################
#### Build pandas dataframes ####
#################################
time            = np.transpose(time)
anomalies       = build_a_bear.anom_to_pandas(
                                                sst_hc,
                                                sst_obs,
                                                time,
                                                lon,
                                                lat,
                                                step0=step0
                                            )

crps_ss         = build_a_bear.score_to_pandas(
                                                crps_ss,
                                                time,
                                                name='CRPS_SS',
                                                step0=step0
                                            )
crps_ss_stat    = build_a_bear.score_to_pandas(
                                                crps_ss_stat,
                                                time,
                                                name='CRPS_SS',
                                                step0=step0
                                            )

acc             = build_a_bear.score_to_pandas(
                                                acc,
                                                time,
                                                name='ACC',
                                                step0=step0
                                            )
# print(
#     anomalies.set_index(
#                         [
#                             'ensemble_member',
#                             'lead_time',
#                             'time',
#                             'month',
#                             'season',
#                             'lon',
#                             'lat'
#                         ]
#                     ).to_numpy().shape
#                 )
# reliability     =
# build_a_bear.pandas_reliability(
#                                                 3,
#                                                 sst_hc,
#                                                 sst_obs,
#                                                 time,
#                                                 lon,
#                                                 lat,
#                                                 step0=step0
#                                             )
# exit()
##################
#### Plotting ####
##################
fr.qq_plot(anomalies[anomalies['lead_time'].isin([10,17,24,31])],domainID)
fr.line_plot(acc,domainID+'_acc',var_name='ACC')
fr.line_plot(crps_ss,domainID+'_crps',var_name='CRPS_SS')
fr.line_plot(crps_ss_stat,domainID+'_crps_stat',var_name='CRPS_SS')
fr.geographic(domainID)
