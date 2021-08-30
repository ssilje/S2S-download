import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd

from .grid import cap_grid
from .handle_domain import get_bounds
from S2S.local_configuration_H import config
import scripts.Henrik.handle_datetime as dt
import scripts.Henrik.figure_receipts as fr

domainID      = 'NVK'
t_start       = (2020,1,23) # start time, fmt: (year,month)
t_end         = (2020,11,15) # end time, fmt: (year,month)

clim_t_start  = (2000,1,1)
clim_t_end    = t_start

from .data_handler import SST

compute_climatology  = 0
compute_fc_anomalies = 0
compute_hc_anomalies = 0

while True:

    if compute_climatology:

        climatology = SST()
        climatology.load('ERA5',clim_t_start,clim_t_end,domainID)

        climpath = '%s%s_%s-%s_%s%s'%(
                                climatology.path['ERA5'],
                                'sst',
                                climatology.time['ERA5'][0].strftime('%Y'),
                                climatology.time['ERA5'][-1].strftime('%Y'),
                                'climatology',
                                '.nc'
                                )

        climatology.ERA5.groupby('time.dayofyear').mean().to_netcdf(climpath)

        compute_climatology = 0

    else:

        climpath = '%s%s_%s-%s_%s%s'%(
                                config['ERA5'],
                                'sst',
                                dt.to_datetime(clim_t_start).strftime('%Y'),
                                dt.to_datetime(clim_t_end).strftime('%Y'),
                                'climatology',
                                '.nc'
                                )

        climatology = xr.open_dataset(climpath)

        break

while True:

    if compute_anomalies:

        sst = SST()

        sst.load('ERA5',t_start,t_end,domainID,match_fc=True)
        sst.load('S2SH',t_start,t_end,domainID)
        sst.load('S2SF',t_start,t_end,domainID)

        obs_anom   = sst.ERA5.groupby('time.dayofyear') - climatology
        model_clim = sst.S2SH.mean('time').mean('number')
        fc_anom    = sst.S2SF - model_clim

        obs_anom = obs_anom.resample(time='7D').mean()
        fc_anom  = fc_anom.resample(step='7D').mean()

        obs_anom = obs_anom.transpose('time','lon','lat')
        fc_anom  = fc_anom.transpose('number','step','cast_time','longitude','latitude')

        fc_anom  = fc_anom.rename(
                                    {
                                        'number':'ensemble_member',
                                        'longitude':'lon',
                                        'latitude':'lat'
                                        }
                                    )

        obs_filepath = '%s%s_%s-%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                'anom_obs',
                                '.nc'
                                )

        fc_filepath  = '%s%s_%s-%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                'fc_obs',
                                '.nc'
                                )

        obs_anom.to_netcdf(obs_filepath)
        fc_anom.to_netcdf(fc_filepath)

        compute_anomalies = 0

    else:

        obs_filepath = '%s%s_%s-%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                'anom_obs',
                                '.nc'
                                )

        fc_filepath  = '%s%s_%s-%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                'fc_obs',
                                '.nc'
                                )

        obs_anom = xr.open_dataset(obs_filepath)
        fc_anom  = xr.open_dataset(fc_filepath)

        break

# Data structures:
#   in:
#       observations    (time,cast_time,lon,lat)
#       forecast        (ensemble,step,cast_time,lon,lat)
#       weights         (ensemble,step,cast_time,lon,lat)
#   out:
#       crps            (time,cast_time,lon,lat)
#
#   validation_time (time) = step + cast_time

validation_time = np.array(
                        [ct+np.array(fc_anom.step)
                            for ct in np.array(fc_anom.time)]
                        )

# add cast_time dimension to observations
observations    = np.stack(
                        [obs_anom.sst.sel(time=t) for t in validation_time],
                        axis=1
                        )

forecast        = np.array(fc_anom.sst)

# weights are given by the cosine of the latitude of each grid point
weights         = np.cos(
                    np.array([
                        np.array([
                            np.array([
                                np.meshgrid(obs_anom.lon,obs_anom.lat)[1]
                                ]*fc_anom.sizes['cast_time'])
                            ]*fc_anom.sizes['step'])
                        ]*fc_anom.sizes['ensemble_member'])
                    )

print(climatology.sst)
exit()
data = SST()
data.load('ERA5',(2020,1,23),(2020,2,15),domainID='NVK',match_fc=True)
data.load('S2SF',(2020,1,23),(2020,2,15),domainID='NVK')
data.load('S2SH',(2020,1,23),(2020,2,15),domainID='NVK')
print(data.ERA5)
print(data.S2SF)
print(data.S2SH)


# fr.qq_plot(forecast,observations)

# print(pd.to_datetime(validation_time[0,0]).dayofyear)
# exit()
# climatology.resample(time='7D').mean()
#
#
# fc_clim   = np.stack(
#                     [climatology.sst.sel(dayofyear=(ct+np.array(fc_anom.step)).dayofyear)
#                         for ct in np.array(fc_anom.time)],axis=1
#                 )
# print(observations.shape,forecast.shape,weights.shape,fc_clim.shape)
#
# crps_fc   = ps.crps_ensemble(
#                 observations=observations,
#                 forecasts=forecast,
#                 weights=weights
#                 )
#
# crps_clim = ps.crps_ensemble(
#                 observations=observations,
#                 forecasts=fc_clim,
#                 weights=weights
#                 )
