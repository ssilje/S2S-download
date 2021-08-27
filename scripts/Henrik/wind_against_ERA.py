import xarray as xr
import pandas as pd
import numpy as np
import xskillscore as xs

import json
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics      import latex, graphics

from matplotlib.colors import BoundaryNorm

from S2S.data_handler  import ERA5, BarentsWatch, Archive
from S2S.process       import Hindcast, Observations, Grid2Point
from S2S               import models
import S2S.scoring     as sc
import S2S.xarray_helpers as xh

t_start = (2019,7,1)
t_end   = (2020,6,25)
clim_t_start  = (1999,1,1)
clim_t_end    = (2021,1,4)
bounds1  = (0,28,35,73)
bounds2  = (345,360,35,73)

steps = pd.to_timedelta([7, 14, 21, 28, 35],'D')
path = '/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/NA/'

# for month in np.arange(1,13,1):
#
#     abs1 = xr.open_dataset(path+'/split_domain/h1_absolute_'+str(month)+'.nc')
#     abs2 = xr.open_dataset(path+'/split_domain/h2_absolute_'+str(month)+'.nc')
#
#     anom1 = xr.open_dataset(path+'/split_domain/h1_anomalies_'+str(month)+'.nc')
#     anom2 = xr.open_dataset(path+'/split_domain/h2_anomalies_'+str(month)+'.nc')
#
#     xr.concat(
#                 [
#                     abs1,
#                     abs2
#                 ],
#             'lon'
#             ).to_netcdf(path+'hindcast_absolute_'+str(month)+'.nc')
#
#     xr.concat(
#                 [
#                     anom1,
#                     anom2
#                 ],
#             'lon'
#             ).to_netcdf(path+'hindcast_anomalies_'+str(month)+'.nc')
#
# exit()
# hindcast_1  = Hindcast(
#                         'U10',
#                         t_start,
#                         t_end,
#                         bounds1,
#                         high_res=False,
#                         steps  =steps,
#                         process=False,
#                         download=False,
#                         split_work=True
#                     )

# hindcast_2  = Hindcast(
#                         'U10',
#                         t_start,
#                         t_end,
#                         bounds2,
#                         high_res=False,
#                         steps  =steps,
#                         process=False,
#                         download=False,
#                         split_work=True
#                     )
#
era1_u = ERA5().load(
                    'u10',
                    clim_t_start,
                    clim_t_end,
                    bounds1,
                    )['u10']

era1_v = ERA5().load(
                    'v10',
                    clim_t_start,
                    clim_t_end,
                    bounds1,
                    )['v10']

era1 = xr.apply_ufunc(
                xh.absolute,era1_u,era1_v,
                input_core_dims  = [[],[]],
                output_core_dims = [[]],
                vectorize=True,dask='parallelized'
            ).rolling(time=7,center=True).mean()

del era1_u
del era1_v

era2_u = ERA5().load(
                    'u10',
                    clim_t_start,
                    clim_t_end,
                    bounds2,
                    )['u10']

era2_v = ERA5().load(
                    'v10',
                    clim_t_start,
                    clim_t_end,
                    bounds2,
                    )['v10']

era2 = xr.apply_ufunc(
                xh.absolute,era2_u,era2_v,
                input_core_dims  = [[],[]],
                output_core_dims = [[]],
                vectorize=True,dask='parallelized'
            ).rolling(time=7,center=True).mean()

del era2_u
del era2_v

xr.concat(
            [
                era1,
                era2
            ],
        'lon'
        ).rename('U10').to_netcdf(path+'era_absolute.nc')


for month in np.arange(1,13,1):

    print(month)

    hc   = xr.open_dataset(path+'hindcast_absolute_'+str(month)+'.nc')['U10']
    time = hc.time
    step = hc.step
    del hc

    obs = xr.open_dataset(path+'era_absolute.nc')['U10']

    print('\tAssign step dimension to observations')
    obs = xh.at_validation(
                                obs,
                                time + step,
                                ddays=1
                                )

    obs = obs.rename('U10')
    obs = obs.drop('validation_time')

    # lon/lat are lost in storage process
    # if contained in encoding['coordinates']
    # coords = obs.encoding['coordinates'].split(' ')
    # while 'lon' in coords: coords.remove('lon')
    # while 'lat' in coords: coords.remove('lat')
    # obs.encoding['coordinates'] = ' '.join(coords)

    obs.to_netcdf(path+'era_absolute_'+str(month)+'.nc')

    print('\tCompute climatology')
    mean,std = xh.o_climatology(obs)

    mean = mean.rename('U10')
    std  = std.rename('U10')

    obs_a = ( obs - mean ) / std

    obs_a.to_netcdf(path+'era_anomalies_'+str(month)+'.nc')
    mean.to_netcdf(path+'era_mean_'+str(month)+'.nc')
    std.to_netcdf(path+'era_std_'+str(month)+'.nc')

    del obs
    del obs_a
    del mean
    del std

    # abs1 = xr.open_dataset(path+'h1_absolute_'+str(month)+'.nc')['U10']
    # abs2 = xr.open_dataset(path+'h2_absolute_'+str(month)+'.nc')['U10']
    #
    # xr.concat(
    #             [
    #                 abs1,
    #                 abs2
    #             ],
    #         'lat'
    #         ).to_netcdf(path+'hindcast_absolute_month_'+str(month)+'.nc')
    #
    # del abs1
    # del abs2
    #
    # anom1 = xr.open_dataset(path+'h1_anomalies_'+str(month)+'.nc')['U10']
    # anom2 = xr.open_dataset(path+'h2_anomalies_'+str(month)+'.nc')['U10']
    #
    # xr.concat(
    #             [
    #                 anom1,
    #                 anom2
    #             ],
    #         'lat'
    #         ).to_netcdf(path+'hindcast_anomalies_month_'+str(month)+'.nc')
    #
    # del anom1
    # del anom2
