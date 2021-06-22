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

from S2S.graphics import mae,crps,brier,spread_error
from S2S import location_cluster

domainID = 'norwegian_coast'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

xlim = 0.5
ylim = 0.25

# observations = BarentsWatch().load('all',no=350).sortby('time')

observations = BarentsWatch().load(['Hisdalen','Stokkvika']).sortby('time')

hindcast     = xh.load_by_location(observations.location,'h_temp')[var]
hindcast_a    = xh.load_by_location(observations.location,'ha_temp')[var]

stacked_obs   = xh.load_by_location(observations.location,'v_temp')[var]
stacked_obs_a = xh.load_by_location(observations.location,'v_a_temp')[var]

clim_mean     = xh.load_by_location(observations.location,'clim_mean_temp')[var]
clim_std      = xh.load_by_location(observations.location,'clim_std_temp')[var]

stacked_era   = xh.load_by_location(observations.location,'era_v_temp')[var]
stacked_era_a = xh.load_by_location(observations.location,'era_v_a_temp')[var]

#
era_mean     = xh.load_by_location(
                                    observations.location,
                                    'era_clim_mean_temp'
                                )[var]
era_std      = xh.load_by_location(
                                    observations.location,
                                    'era_clim_std_temp'
                                )[var]

# hindcast      = hindcast.sel(
#                             time=slice(stacked_obs.time.min(),
#                             stacked_obs.time.max())
#                             )

clim_fc       = models.clim_fc(clim_mean,clim_std)
era_fc        = models.clim_fc(era_mean,era_std)

hindcast      = xh.assign_validation_time(hindcast)
hindcast_a    = xh.assign_validation_time(hindcast_a)

clim_mean     = xh.assign_validation_time(clim_mean)
clim_std      = xh.assign_validation_time(clim_std)
clim_fc       = xh.assign_validation_time(clim_fc)
stacked_obs   = xh.assign_validation_time(stacked_obs)
stacked_obs_a = xh.assign_validation_time(stacked_obs_a)

era_mean      = xh.assign_validation_time(era_mean)
era_std       = xh.assign_validation_time(era_std)
era_fc        = xh.assign_validation_time(era_fc)
stacked_era   = xh.assign_validation_time(stacked_era)
stacked_era_a = xh.assign_validation_time(stacked_era_a)


hindcast_ba = models.bias_adjustment_torrabla(
                                            hindcast_a,
                                            stacked_era_a
                                            )
gr.timeseries(
                stacked_era,
                cast=[hindcast_a*era_std + era_mean,hindcast_ba+era_mean],
                title='EC ba',
                filename='era_ba',
                clabs=['EC','EC ba']
            )
#
crps.skill_agg(
            stacked_era_a,
            [hindcast_a,hindcast_ba],
            xr.full_like(stacked_era_a,0),
            xr.full_like(stacked_era_a,1),
            title='ERA Simple Bias Adjustment',
            filename='crps_ERA_ba',
            dim='validation_time.month',
            mlabs=['sb','ba']
        )

hindcast_ba = models.bias_adjustment_torrabla(
                                            hindcast_a,
                                            stacked_obs_a
                                            )
gr.timeseries(
                stacked_obs_a,
                cast=[hindcast_a,hindcast_ba],
                title='EC ba',
                filename='BW_ba',
                clabs=['EC','EC ba']
            )
#
crps.skill_agg(
            stacked_obs_a,
            [hindcast_a,hindcast_ba],
            xr.full_like(stacked_era_a,0),
            xr.full_like(stacked_era_a,1),
            title='ERA Bias Adjustment',
            filename='crps_BW_ba',
            dim='validation_time.month',
            mlabs=['sb','ba']
        )
exit()

# gr.point_map(hindcast,poi='Hisdalen')

# gr.timeseries(
#                 stacked_era,
#                 cast=[era_fc,hindcast],
#                 title='EC',
#                 filename='era',
#                 clabs=['clim','EC']
#             )
#
# gr.timeseries(
#                 stacked_era,
#                 cast=[era_fc,hindcast_a*era_std + era_mean],
#                 title='EC',
#                 filename='era_sb',
#                 clabs=['clim','EC']
#             )
#
# gr.timeseries(
#                 stacked_era,
#                 cast=[era_fc,hindcast],
#                 title='EC',
#                 filename='era',
#                 clabs=['clim','EC'],
#                 lead_time=[30,37]
#             )
#
# gr.timeseries(
#                 stacked_era,
#                 cast=[era_fc,hindcast_a*era_std + era_mean],
#                 title='EC',
#                 filename='era_sb',
#                 clabs=['clim','EC'],
#                 lead_time=[30,37]
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,hindcast],
#                 title='EC',
#                 filename='barentswatch',
#                 clabs=['clim','EC']
#             )
#
# gr.timeseries(
#                 stacked_obs_a,
#                 cast=[hindcast_a],
#                 title='EC',
#                 filename='barentswatch_a',
#                 clabs=['EC']
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,(hindcast_a*clim_std)+clim_mean],
#                 title='EC',
#                 filename='barentswatch_sb',
#                 clabs=['clim','EC']
#              )
#
# gr.timeseries(
#                 stacked_obs_a,
#                 cast=[hindcast_a],
#                 title='EC',
#                 filename='barentswatch_a',
#                 clabs=['EC']
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,(hindcast_a*clim_std)+clim_mean],
#                 title='EC',
#                 filename='barentswatch_sb',
#                 clabs=['clim','EC']
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,hindcast],
#                 title='EC',
#                 filename='barentswatch',
#                 clabs=['clim','EC'],
#                 lead_time = [30,37]
#             )
#
# gr.timeseries(
#                 stacked_obs_a,
#                 cast=[hindcast_a],
#                 title='EC',
#                 filename='barentswatch_a',
#                 clabs=['EC'],
#                 lead_time=[30,37]
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,(hindcast_a*clim_std)+clim_mean],
#                 title='EC',
#                 filename='barentswatch_sb',
#                 clabs=['clim','EC'],
#                 lead_time=[30,37]
#              )
#
# gr.timeseries(
#                 stacked_obs_a,
#                 cast=[hindcast_a],
#                 title='EC',
#                 filename='barentswatch_a',
#                 clabs=['EC'],
#                 lead_time=[30,37]
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,(hindcast_a*clim_std)+clim_mean],
#                 title='EC',
#                 filename='barentswatch_sb',
#                 clabs=['clim','EC'],
#                 lead_time=[30,37]
#             )
#
mae.skill_agg(
            stacked_era,
            [hindcast],
            era_mean,
            era_std,
            title='ERA',
            filename='mae_ERA',
            dim='validation_time.month'
        )

mae.skill_agg(
            stacked_era,
            [(hindcast_a*era_std)+era_mean],
            era_mean,
            era_std,
            title='ERA Simple Bias Adjustment',
            filename='mae_ERA_sb',
            dim='validation_time.month'
        )
mae.skill_agg(
            stacked_era_a,
            [hindcast_a],
            xr.full_like(stacked_era_a,0),
            xr.full_like(stacked_era_a,1),
            title='ERA anomalies',
            filename='mae_ERA_a',
            dim='validation_time.month'
        )

mae.skill_agg(
            stacked_obs,
            [hindcast],
            clim_mean,
            clim_std,
            title='Barentswatch',
            filename='mae_BW',
            dim='validation_time.month'
        )
#
mae.skill_agg(
            stacked_obs,
            [(hindcast_a*clim_std)+clim_mean],
            clim_mean,
            clim_std,
            title='Barentswatch Simple Bias Adjustment',
            filename='mae_BW_sb',
            dim='validation_time.month'
        )

mae.skill_agg(
            stacked_obs_a,
            [hindcast_a],
            xr.full_like(stacked_obs_a,0),
            xr.full_like(stacked_obs_a,1),
            title='Barentswatch anomalies',
            filename='mae_BW_a',
            dim='validation_time.month'
        )

crps.skill_agg(
            stacked_era,
            [hindcast],
            era_mean,
            era_std,
            title='ERA',
            filename='crps_ERA',
            dim='validation_time.month'
        )

crps.skill_agg(
            stacked_era,
            [(hindcast_a*era_std)+era_mean],
            era_mean,
            era_std,
            title='ERA Simple Bias Adjustment',
            filename='crps_ERA_sb',
            dim='validation_time.month'
        )

crps.skill_agg(
            stacked_era_a,
            [hindcast_a],
            xr.full_like(stacked_era_a,0),
            xr.full_like(stacked_era_a,1),
            title='ERA anomalies',
            filename='crps_ERA_a',
            dim='validation_time.month'
        )

crps.skill_agg(
            stacked_obs,
            [hindcast],
            clim_mean,
            clim_std,
            title='Barentswatch',
            filename='crps_BW',
            dim='validation_time.month'
        )

crps.skill_agg(
            stacked_obs,
            [(hindcast_a*clim_std)+clim_mean],
            clim_mean,
            clim_std,
            title='Barentswatch Simple Bias Adjustment',
            filename='crps_BW_sb',
            dim='validation_time.month'
        )
crps.skill_agg(
            stacked_obs_a,
            [hindcast_a],
            xr.full_like(stacked_obs_a,0),
            xr.full_like(stacked_obs_a,1),
            title='Barentswatch anomalies',
            filename='crps_BW_a',
            dim='validation_time.month'
        )
#
# brier.skill_agg(
#             stacked_era,
#             [hindcast_a*era_std + era_mean],
#             era_mean,
#             era_std,
#             std=1.5,
#             dim='validation_time.month',
#             filename='ERA_brier_sb',
#             title='ERA Simple Bias Adjustment',
#             ylab='BrierSS'
#         )
#
# brier.skill_agg(
#             stacked_obs,
#             [hindcast_a*clim_std + clim_mean],
#             clim_mean,
#             clim_std,
#             std=2.,
#             dim='validation_time.month',
#             filename='BW_brier_sb',
#             title='Barentswatch Simple Bias Adjustment',
#             ylab='BrierSS'
#         )

# mae.map(
#         stacked_era,
#         hindcast_a*era_std + era_mean,
#         era_mean,
#         era_std,
#         dim='validation_time.month',
#         title='ERA Simple Bias Adjustment',
#         filename='ERA_sb'
#         )
#
# mae.map(
#         stacked_obs,
#         hindcast_a*clim_std + clim_mean,
#         clim_mean,
#         clim_std,
#         dim='validation_time.season',
#         title='Barentswatch Simple Bias Adjustment',
#         filename='BW_sb'
#         )

# spread_error.r_ss(
#             stacked_era,
#             [hindcast,hindcast_a*era_std + era_mean],
#             era_mean,
#             era_std,
#             labels=['EC','EC sb'],
#             dim='validation_time.month',
#             filename='ERA_se_sb',
#             title='ERA Simple Bias Adjustment',
#             ylab='Spread Error SS'
#         )
#
# spread_error.r_ss(
#             stacked_obs,
#             [hindcast,hindcast_a*clim_std + clim_mean],
#             clim_mean,
#             clim_std,
#             labels=['EC','EC sb'],
#             dim='validation_time.month',
#             filename='BW_se_sb',
#             title='Barentswatch Simple Bias Adjustment',
#             ylab='Spread Error SS'
#         )

# mae.cluster_map(
#         stacked_era,
#         hindcast_a*era_std + era_mean,
#         era_mean,
#         era_std,
#         c_lim=(xlim,ylim),
#         dim='validation_time.month',
#         title='ERA Simple Bias Adjustment',
#         filename='ERA_sb'
#         )
#
# mae.cluster_map(
#         stacked_obs,
#         hindcast_a*clim_std + clim_mean,
#         clim_mean,
#         clim_std,
#         c_lim=(xlim,ylim),
#         dim='validation_time.month',
#         title='Barentswatch Simple Bias Adjustment',
#         filename='BW_sb'
#         )
#
# exit()
# xlim = 0.5
# ylim = 0.25
#
# hindcast   = location_cluster.cluster(hindcast,'Hisdalen',xlim,ylim)
# hindcast_a = location_cluster.cluster(hindcast_a,'Hisdalen',xlim,ylim)
#
# clim_fc       = location_cluster.cluster(clim_fc,'Hisdalen',xlim,ylim)
# stacked_obs   = location_cluster.cluster(stacked_obs,'Hisdalen',xlim,ylim)
# stacked_obs_a = location_cluster.cluster(stacked_obs_a,'Hisdalen',xlim,ylim)
#
# era_fc        = location_cluster.cluster(era_fc,'Hisdalen',xlim,ylim)
# stacked_era   = location_cluster.cluster(stacked_era,'Hisdalen',xlim,ylim)
# stacked_era_a = location_cluster.cluster(stacked_era_a,'Hisdalen',xlim,ylim)
#
# gr.point_map(hindcast,poi='Hisdalen')
#
# mae.skill_agg_cluster(
#             stacked_era,
#             [hindcast],
#             era_mean,
#             era_std,
#             loc='Hisdalen',
#             title='ERA',
#             dim='validation_time.season',
#             filename='mae_cluster_ERA'
#         )
#
# mae.skill_agg_cluster(
#             stacked_era,
#             [hindcast_a*era_std + era_mean],
#             era_mean,
#             era_std,
#             loc='Hisdalen',
#             title='ERA Simple Bias Adjustment',
#             dim='validation_time.season',
#             filename='mae_cluster_ERA_sb'
#         )
#
# mae.skill_agg_cluster(
#             stacked_obs,
#             [hindcast],
#             clim_mean,
#             clim_std,
#             loc='Hisdalen',
#             title='Barentswatch',
#             dim='validation_time.season',
#             filename='mae_cluster_BW'
#         )
#
# mae.skill_agg_cluster(
#             stacked_obs,
#             [hindcast_a*clim_std + clim_mean],
#             clim_mean,
#             clim_std,
#             loc='Hisdalen',
#             title='Barentswatch Simple Bias Adjustment',
#             dim='validation_time.season',
#             filename='mae_cluster_BW_sb'
#         )
# exit()
################################################################################
# mae.plot(
#             stacked_obs,
#             [(hindcast_a*clim_std)+clim_mean],
#             clim_mean,
#             clim_std
#         )
# exit()
# crps_fc   = sc.crps_ensemble(stacked_obs,hindcast)
# crps_clim = xs.crps_gaussian(
#                         stacked_obs,
#                         clim_mean,
#                         clim_std.where(clim_std!=0),
#                         dim=[]
#                         )
#
# rmse_fc   = xs.rmse(stacked_obs,hindcast.mean('member'),dim=[])
# rmse_clim = xs.rmse(stacked_obs,clim_mean.reindex(time=stacked_obs.time),dim=[])
#
# # gr.residual_plot(hindcast,stacked_obs)
# # exit()
# gr.skill_plot(
#                 rmse_fc,
#                 rmse_clim,
#                 title='EC',
#                 filename='RMSEss_EC',
#                 ylab='RMSESS',
#                 dim='validation_time.month'
#             )
#
# gr.skill_plot(
#                 crps_fc,
#                 crps_clim,
#                 title='EC',
#                 filename='CRPss_EC',
#                 ylab='CRPSS',
#                 dim='validation_time.month'
#             )
# ############################################################################
# crps_fc   = sc.crps_ensemble(stacked_obs_a,hindcast_a)
# crps_clim = xs.crps_gaussian(
#                         stacked_obs_a,
#                         0,
#                         1,
#                         dim=[]
#                         )
#
# rmse_fc   = xs.rmse(stacked_obs_a,hindcast_a.mean('member'),dim=[])
# rmse_clim = xs.rmse(stacked_obs_a,xr.full_like(stacked_obs_a,0),dim=[])
#
# gr.skill_plot(
#                 rmse_fc,
#                 rmse_clim,
#                 title='EC simple bias adjustment',
#                 filename='RMSEss_EC_a',
#                 ylab='RMSESS',
#                 dim='validation_time.month'
#             )
#
# gr.skill_plot(
#                 crps_fc,
#                 crps_clim,
#                 title='EC simple bias adjustment',
#                 filename='CRPss_EC_a',
#                 ylab='CRPSS',
#                 dim='validation_time.month'
#             )
# ############################################################################
# crps_fc   = sc.crps_ensemble(stacked_obs,(hindcast_a*clim_std)+clim_mean)
# crps_clim = xs.crps_gaussian(
#                         stacked_obs,
#                         clim_mean,
#                         clim_std,
#                         dim=[]
#                         )
#
# rmse_fc   = xs.rmse(stacked_obs,((hindcast_a*clim_std)+clim_mean).mean('member'),dim=[])
# rmse_clim = xs.rmse(stacked_obs,clim_mean,dim=[])
#
# gr.skill_plot(
#                 rmse_fc,
#                 rmse_clim,
#                 title='EC simple bias adjustment',
#                 filename='RMSEss_EC_simple',
#                 ylab='RMSESS',
#                 dim='validation_time.month'
#             )
#
# gr.skill_plot(
#                 crps_fc,
#                 crps_clim,
#                 title='EC simple bias adjustment',
#                 filename='CRPss_EC_simple',
#                 ylab='CRPSS',
#                 dim='validation_time.month'
#             )
# ################################################################################
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,hindcast],
#                 title='EC',
#                 filename='test',
#                 clabs=['clim','EC']
#             )
#
# gr.timeseries(
#                 stacked_obs_a,
#                 cast=[hindcast_a],
#                 title='EC',
#                 filename='test_anomaly',
#                 clabs=['EC']
#             )
#
# gr.timeseries(
#                 stacked_obs,
#                 cast=[clim_fc,(hindcast_a*clim_std)+clim_mean],
#                 title='EC',
#                 filename='test_simple',
#                 clabs=['clim','EC']
#             )
# exit()
