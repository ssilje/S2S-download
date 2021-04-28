import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd

from .handle_domain import get_bounds
from S2S.local_configuration import config
import scripts.Henrik.handle_datetime as dt
import scripts.Henrik.figure_receipts as fr
from .acc import ACC
from .create_domain_file import make_dir
import scripts.Henrik.build_a_bear as build_a_bear

domainID      = 'NVS'

t_start  = (2020,1,23)
t_end    = (2021,1,14)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

from .data_handler import SST

get_observations = 0
get_hindcast = 0

while True:

    if get_observations:

        observations = SST()
        observations.load('ERA5',clim_t_start,clim_t_end,domainID)

        obspath = '%s%s_%s-%s_%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(clim_t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(clim_t_end).strftime('%Y-%m-%d'),
                                'reanalysis',
                                domainID,
                                '.nc'
                                )
        make_dir('/'.join(obspath.split('/')[:-1]))
        observations.ERA5.transpose('time','lon','lat').to_netcdf(obspath)

        get_observations = 0

    else:

        obspath = '%s%s_%s-%s_%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(clim_t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(clim_t_end).strftime('%Y-%m-%d'),
                                'reanalysis',
                                domainID,
                                '.nc'
                                )

        observations = xr.open_dataset(obspath)

        break

while True:

    if get_hindcast:

        sst = SST()
        sst.load('S2SH',t_start,t_end,domainID)

        hindcast = sst.S2SH.transpose(
                                    'number','step','time',
                                    'longitude','latitude'
                                    )

        hindcast = hindcast.rename(
                                    {
                                        'number':'ensemble_member',
                                        'longitude':'lon',
                                        'latitude':'lat'
                                        }
                                    )

        hc_filepath = '%s%s_%s-%s_%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_end).strftime('%Y-%m-%d'),
                                'hc',
                                domainID,
                                '.nc'
                                )

        make_dir('/'.join(hc_filepath.split('/')[:-1]))
        hindcast.to_netcdf(hc_filepath)

        get_hindcast = 0

    else:

        hc_filepath = '%s%s_%s-%s_%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_end).strftime('%Y-%m-%d'),
                                'hc',
                                domainID,
                                '.nc'
                                )

        hindcast = xr.open_dataset(hc_filepath)

        break

####################################################
#### Cap time of hindcast to match observations ####
####################################################
time_index = np.argmin(
                np.abs(
                    np.array(
                        hindcast.time+hindcast.step[-1]-observations.time[-1]
                        )
                    )
                ) - 1
hindcast = hindcast.isel(time=slice(0,time_index))

#####################################
#### Compute anomalies and modes ####
#####################################
obs_mean = observations.groupby('time.dayofyear').mean() # not necissary to store this variable
obs_std  = observations.groupby('time.dayofyear').std()
obs_anom = observations.groupby('time.dayofyear') - obs_mean

hc_mean  = hindcast.groupby('time.dayofyear').mean().mean('ensemble_member') # not necissary to store this variable
hc_anom  = hindcast.groupby('time.dayofyear') - hc_mean

############################################
#### Compute running means of 7D window ####
############################################
obs_std  = obs_std.sst\
            .pad(dayofyear=7,mode='wrap')\
                .rolling(dayofyear=7,center=True).mean()\
                    .isel(dayofyear=slice(7,-7))\
                        .to_dataset(name='sst')

obs_anom = obs_anom.sst\
            .rolling(time=7,center=True).mean()\
                .to_dataset(name='sst')

hc_anom  = hc_anom.sst\
            .rolling(step=7,center=True).mean()\
                .dropna('step')\
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

#############################
#### Compute CRPS and SS ####
#############################
crps_fc   = ps.crps_ensemble(
                observations=sst_obs,
                forecasts=np.moveaxis(sst_hc,0,-1),
                weights=np.moveaxis(weights,0,-1)
                )

crps_clim = ps.crps_gaussian(x=sst_obs, mu=0, sig=std_obs) # does not weight with cos(lat), might be problematic in next operation crps_SS

crps_ss   = 1 - np.nanmean(crps_fc,axis=(-1,-2))/np.nanmean(crps_clim,axis=(-1,-2))

#############
#### ACC ####
#############
acc = ACC(sst_hc,sst_obs,weights)

#################################
#### Build pandas dataframes ####
#################################
time      = np.transpose(time)
anomalies = build_a_bear.anom_to_pandas(sst_hc,sst_obs,time,lon,lat,step0=step0)
crps_ss   = build_a_bear.score_to_pandas(crps_ss,time,name='CRPS_SS',step0=step0)
acc       = build_a_bear.score_to_pandas(acc,time,name='ACC',step0=step0)

##################
#### Plotting ####
##################
fr.qq_plot(anomalies[anomalies['lead_time'].isin([10,17,24,31])],domainID)
fr.line_plot(acc,domainID+'_acc',var_name='ACC')
fr.line_plot(crps_ss,domainID+'_crps',var_name='CRPS_SS')
fr.geographic(domainID)
