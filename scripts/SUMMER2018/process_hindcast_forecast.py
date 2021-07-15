import pandas as pd
import xarray as xr
import xskillscore as xs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from S2S.data_handler import ERA5
from S2S.process import Hindcast, Observations, Forecast

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models
import S2S.xarray_helpers    as xh
from S2S.scoring import uncentered_acc, centered_acc
from S2S.local_configuration import config




bounds          = (0,28,55,75)

var             = 't2m'
clabel          = 'K'
#var             = 'abs_wind'
t_start         = (2018,3,1)
t_end           = (2018,9,24)


clim_t_start    = (1999,1,1)
clim_t_end      = (2019,12,31)


high_res        = False

steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

print('\tProcess Forecast')
grid_hindcast = Forecast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        steps=steps,
                        process=True,
                        download=False,
                        split_work=True
                    )


print('\tProcess Hindcast')
grid_hindcast = Hindcast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        steps=steps,
                        process=True,
                        download=False,
                        split_work=True
                    )


dim             = 'validation_time.month'



print('\tLoad ERA')
if var == 'abs_wind':
    era_u = ERA5(high_res=high_res)\
                         .load('u10',clim_t_start,clim_t_end,bounds)['u10']

    era_v = ERA5(high_res=high_res)\
                        .load('v10',clim_t_start,clim_t_end,bounds)['v10']

    era_u,era_v = xr.align(era_u,era_v)

    era = xr.apply_ufunc(
                    xh.absolute,era_u,era_v,
                    input_core_dims  = [[],[]],
                    output_core_dims = [[]],
                    vectorize=True,dask='parallelized'
                )

    era = era.rename(var)

else:

    era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]



grid_observations = Observations(
                            name='Era',
                            observations=era,
                            forecast=grid_hindcast,
                            process=False
                            )

