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

observations = BarentsWatch().load('all',no=350).sortby('time')

# val_obs      = xh.load_by_location(observations.location,'v_temp')
# val_obs_a    = xh.load_by_location(observations.location,'v_a_temp')
clim_mean    = xh.load_by_location(observations.location,'clim_mean_temp')
clim_std     = xh.load_by_location(observations.location,'clim_std_temp')
# hindcast_a     = xh.load_by_location(observations.location,'hp_temp')

def cluster(da,loc,lon_tolerance,lat_tolerance):

    z = da.sel(location=loc)

    lon_lim = np.array([-lon_tolerance,lon_tolerance]) + z.lon.values
    lat_lim = np.array([-lat_tolerance,lat_tolerance]) + z.lat.values

    c = da.where(lon_lim[0]<da.lon,drop=True)
    c = c.where(c.lon<lon_lim[1],drop=True)

    c = c.where(lat_lim[0]<c.lat,drop=True)
    c = c.where(c.lat<lat_lim[1],drop=True)

    return c

def cluster_dist(
                    hindcast,
                    observations,
                    clim_m,
                    dim='validation_time.month',
                    dist_func=xs.rmse
                ):
    print('\tcluster_mae()')

    N   = len(hindcast.location)
    fc   = []
    clim = []
    for n,loc in enumerate(hindcast.location):

        xh.print_progress(n,N)

        hc = cluster(hindcast,loc,0.5,0.25)
        oc = cluster(observations,loc,0.5,0.25)

        fc_score   = dist_func(oc,hc.mean('member',skipna=True),dim=[])
        clim_score = dist_func(oc,clim_mean,dim=[])

        fc_score = fc_score.mean('location',skipna=True).assign_coords(
                        {
                            'location':loc,
                            'lon':hindcast.lon.sel(location=loc),
                            'lat':hindcast.lat.sel(location=loc)
                        }
                    )
        clim_score = clim_score.mean('location',skipna=True).assign_coords(
                        {
                            'location':loc,
                            'lon':hindcast.lon.sel(location=loc),
                            'lat':hindcast.lat.sel(location=loc)
                        }
                    )

        fc.append(fc_score)
        clim.append(clim_score)

    return xr.concat(fc,'location'),xr.concat(clim,'location')

def cluster_crps(
                    hindcast,
                    observations,
                    clim_m,
                    clim_s,
                    dim='validation_time.month'
                ):
    print('\tcluster_mae()')

    N   = len(hindcast.location)
    fc   = []
    clim = []
    for n,loc in enumerate(hindcast.location):

        xh.print_progress(n,N)

        hc = cluster(hindcast,loc,0.5,0.25)
        oc = cluster(observations,loc,0.5,0.25)

        fc_score   = xs.crps_ensemble(oc,hc,dim=[])
        clim_score = xs.crps_(oc,xr.full_like(oc,0),dim=[])

        fc_score = fc_score.mean('location',skipna=True).assign_coords(
                        {
                            'location':loc,
                            'lon':hindcast.lon.sel(location=loc),
                            'lat':hindcast.lat.sel(location=loc)
                        }
                    )
        clim_score = clim_score.mean('location',skipna=True).assign_coords(
                        {
                            'location':loc,
                            'lon':hindcast.lon.sel(location=loc),
                            'lat':hindcast.lat.sel(location=loc)
                        }
                    )

        fc.append(fc_score)
        clim.append(clim_score)

    return xr.concat(fc,'location'),xr.concat(clim,'location')


# observations = BarentsWatch().load('all',no=350).sortby('time')

val_obs      = xh.load_by_location(observations.location,'v_temp')
hindcast_abs = xh.load_by_location(observations.location,'h_temp')
# val_obs_a    = xh.load_by_location(observations.location,'v_a_temp')
hindcast     = xh.load_by_location(observations.location,'hp_temp')

hindcast_abs = hindcast_abs.sel(time=slice(val_obs.time.min(),val_obs.time.max()))
hindcast     = hindcast.sel(time=slice(val_obs.time.min(),val_obs.time.max()))

rmse_fc,rmse_clim = cluster_dist(hindcast*clim_std + clim_mean,val_obs,clim_mean)
armse_fc,armse_clim = cluster_dist(hindcast_abs,val_obs,clim_mean)

xh.store_by_location(rmse_fc,'rmse_fc_temp')
xh.store_by_location(rmse_clim,'rmse_clim_temp')

xh.store_by_location(armse_fc,'armse_fc_temp')
xh.store_by_location(armse_clim,'armse_clim_temp')

rmse_fc = xh.load_by_location(observations.location,'rmse_fc_temp')
rmse_clim = xh.load_by_location(observations.location,'rmse_clim_temp')

armse_fc = xh.load_by_location(observations.location,'armse_fc_temp')
armse_clim = xh.load_by_location(observations.location,'armse_clim_temp')

# gr.skill_map(
#             rmse_fc[var],
#             rmse_clim[var],
#             title='RMSE SS',
#             filename='RMSEss_map',
#             lead_time=[7,14,21,28,35,42]
#             )

gr.skill_map(
            rmse_fc[var],
            rmse_clim[var],
            title='RMSE SS',
            filename='RMSEss_map',
            lead_time=[7,14,21,28,35,42]
            )

gr.skill_plot(
                rmse_fc[var],
                rmse_clim[var],
                title='EC',
                filename='RMSESS_EC',
                ylab='RMSESS',
                dim='validation_time.month'
            )

gr.skill_map(
            armse_fc[var],
            armse_clim[var],
            title='RMSE SS',
            filename='RMSEss_map',
            lead_time=[7,14,21,28,35,42]
            )

gr.skill_plot(
                armse_fc[var],
                armse_clim[var],
                title='EC',
                filename='RMSESS_EC',
                ylab='RMSESS',
                dim='validation_time.month'
            )

# # probabilistic skill
# crps_mod  = xs.crps_ensemble(val_obs_a[var],hindcast[var],dim=[])
# crps_clim = xs.crps_gaussian(
#                             val_obs_a[var],
#                             mu=0,
#                             sig=1,
#                             dim=[]
#                             )
#
# gr.skill_plot(
#                 crps_mod,
#                 crps_clim,
#                 title='EC',
#                 filename='CRPSS_EC',
#                 ylab='CRPSS'
#             )

# # deterministic skill
# mae_mod  = xs.mae(val_obs_a[var],hindcast[var].mean('member'),dim=[])
# mae_clim = xs.mae(
#                   val_obs_a[var],
#                   xr.full_like(hindcast[var].mean('member'),0),
#                   dim=[]
#                  )
#
# gr.skill_plot(
#                 mae_mod,
#                 mae_clim,
#                 title='EC',
#                 filename='MAESS_EC',
#                 ylab='MAESS',
#                 dim='validation_time.year'
#             )

# rmse_mod  = xs.rmse(val_obs_a[var],hindcast[var].mean('member'),dim=[])
# rmse_clim = xs.rmse(
#                   val_obs_a[var],
#                   xr.full_like(hindcast[var].mean('member'),0),
#                   dim=[]
#                  )
#
# gr.skill_plot(
#                 mae_mod,
#                 mae_clim,
#                 title='EC',
#                 filename='RMSESS_EC',
#                 ylab='RMSESS'
#             )
#
# gr.qq_plot(val_obs_a[var],hindcast[var].mean('member'),
#                             y_axlabel='EC',filename='EC')
#
# hindcast   = (hindcast*clim_std) + clim_mean
# clim_fc    = models.clim_fc(clim_mean,clim_std)
#
#
# gr.timeseries(
#                 val_obs[var],
#                 cast=[clim_fc[var],hindcast[var]],
#                 title='EC',
#                 filename='EC',
#                 clabs=['clim','EC']
#             )
# gr.timeseries(
#                 val_obs[var],
#                 cast=[clim_fc[var],hindcast[var]],
#                 lead_time=[28,35,42],
#                 title='EC',
#                 filename='EC',
#                 clabs=['clim','EC']
#             )
