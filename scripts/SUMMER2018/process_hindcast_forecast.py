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
grid_forecast = Forecast(
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
                        process=False,
                        download=False,
                        split_work=True
                    )

era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]
grid_observations = Observations(
                            name='Era',
                            observations=era,
                            forecast=grid_hindcast,
                            process=False
                            )



stacked_era       = xh.assign_validation_time(grid_observations.data)
stacked_era_a     = xh.assign_validation_time(grid_observations.data_a)
stacked_era_mean  = xh.assign_validation_time(grid_observations.mean)
stacked_era_std   = xh.assign_validation_time(grid_observations.std)

hindcast          = xh.assign_validation_time(grid_hindcast.data)
hindcast_a        = xh.assign_validation_time(grid_hindcast.data_a)
hindcast_mean     = xh.assign_validation_time(grid_hindcast.mean)
hindcast_std      = xh.assign_validation_time(grid_hindcast.std)

forecast          = xh.assign_validation_time(grid_forecast.data)
forecast_a        = forecast - hindcast_mean
#forecast_a        = xh.assign_validation_time(grid_forecast.data_a)
#forecast_mean     = xh.assign_validation_time(grid_forecast.mean)
#forecast_std      = xh.assign_validation_time(grid_forecast.std)



print('\tLoop through all steps')

for lt in steps:

    mod         = model.sel(step=pd.Timedelta(lt,'D'))
    obs         = observations.sel(step=pd.Timedelta(lt,'D'))
    cm          = xr.full_like(observations,0).sel(step=pd.Timedelta(lt,'D')) ## er null pga bruker anomalier
    obs_random  = random_fc_a.sel(step=pd.Timedelta(lt,'D'))

    era         = stacked_era.sel(step=pd.Timedelta(lt,'D'))
    hc          = hindcast.sel(step=pd.Timedelta(lt,'D'))
    
   
  
    x_group     = list(mod.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    y_group     = list(obs.groupby(dim))
    cm_group    = list(cm.groupby(dim))
    yr_group    = list(obs_random.groupby(dim))

    era_group   = list(era.groupby(dim))
    hc_group    = list(hc.groupby(dim))
  
    mae         = []
    c           = [] #lagar en ny xarray med score for kvar mnd
    acc         = [] #lagar en ny xarray med ACC for kvar mnd
    
    era_tmp     = []
    hc_tmp      = []    





