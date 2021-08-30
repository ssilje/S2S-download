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
import scripts.Henrik.compute_reliability as cr
from scripts.Henrik.prepare_data import get_observations, get_hindcast

domainID = 'NVK'

t_start  = (2020,1,23)
t_end    = (2021,1,14)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

###################
#### Load data ####
###################
# observations = get_observations(domainID,clim_t_start,clim_t_end,download=0)
hindcast = get_hindcast(domainID,t_start,t_end)
exit()

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
