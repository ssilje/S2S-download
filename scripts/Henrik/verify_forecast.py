import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd

from .grid import cap_grid
from .handle_domain import get_bounds
from S2S.local_configuration_H import config
import scripts.Henrik.handle_datetime as dt
import scripts.Henrik.figure_receipts as fr
from .acc import ACC
from .create_domain_file import make_dir
import scripts.Henrik.build_a_bear as build_a_bear

domainID      = 'NVS'

t_start  = (2020,1,23)
t_end    = (2020,11,15)

clim_t_start  = (2000,1,23)
clim_t_end    = t_end

from .data_handler import SST

get_observations  = 0
compute_anomalies = 0

while True:

    if get_observations:

        observations = SST()
        observations.load('ERA5',clim_t_start,clim_t_end,domainID,match_fc=True)

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

    if compute_anomalies:

        sst = SST()
        sst.load('S2SH',t_start,t_end,domainID)

        hc_anom   = sst.S2SH.groupby('time.dayofyear') -\
                        sst.S2SH.groupby('time.dayofyear').mean().mean('number')

        hc_anom = hc_anom.resample(step='7D').mean()

        hc_anom = hc_anom.transpose('number','step','time','longitude','latitude')

        hc_anom = hc_anom.rename(
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
                                'anom_hc',
                                domainID,
                                '.nc'
                                )

        make_dir('/'.join(hc_filepath.split('/')[:-1]))
        hc_anom.to_netcdf(hc_filepath)

        compute_anomalies = 0

    else:

        hc_filepath = '%s%s_%s-%s_%s_%s%s'%(
                                config['VALID_DB'],
                                'sst',
                                dt.to_datetime(t_start).strftime('%Y-%m-%d'),
                                dt.to_datetime(t_end).strftime('%Y-%m-%d'),
                                'anom_hc',
                                domainID,
                                '.nc'
                                )

        hc_anom = xr.open_dataset(hc_filepath)

        break


lon = np.array(hc_anom.lon)
lat = np.array(hc_anom.lat)
validation_time = np.array(
                        [ct+np.array(hc_anom.step)
                            for ct in np.array(hc_anom.time)]
                        )

obs_anom        = observations.groupby('time.dayofyear') -\
                    observations.groupby('time.dayofyear').mean()

# stack to match validation time
sst_obs         = np.stack(
                        [obs_anom.sst.sel(time=t) for t in validation_time],
                        axis=1
                        )

sst_hc          = np.array(hc_anom.sst)

# weights are given by the cosine of the latitude of each grid point
weights         = np.cos(
                    np.deg2rad(
                        np.array([
                            np.array([
                                np.array([
                                    np.meshgrid(obs_anom.lat,obs_anom.lon)[0]
                                ]*hc_anom.sizes['time'])
                            ]*hc_anom.sizes['step'])
                        ]*hc_anom.sizes['ensemble_member'])
                    )
                )

# continous ranked probability score
crps_fc   = ps.crps_ensemble(
                observations=sst_obs,
                forecasts=np.moveaxis(sst_hc,0,-1),
                weights=np.moveaxis(weights,0,-1)
                )

# anomaly correlation coefficients
acc = ACC(sst_hc,sst_obs,weights)

# convert to pandas dataframe
time = np.transpose(validation_time)

anomalies = build_a_bear.anom_to_pandas(sst_hc,sst_obs,time,lon,lat)
crps      = build_a_bear.score_to_pandas(crps_fc,time,lon,lat)
acc       = build_a_bear.score_to_pandas(acc,time,name='ACC')

fr.qq_plot(anomalies[anomalies['lead_time']<=4],domainID)
fr.line_plot(acc[acc['lead_time']<=4],domainID,var_name='ACC')
# quantile plot
# fr.qq_plot(sst_hc,sst_obs,domainID)
