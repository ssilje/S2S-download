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
t_start         = (2018,3,1)
t_end           = (2018,9,24)


clim_t_start    = (1999,1,1)
clim_t_end      = (2019,12,31)


high_res        = False

#steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')
steps           = pd.to_timedelta(np.linspace(1,46,46),'D')

print('\tProcess Forecast')
grid_forecast = Forecast(
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
hc_fc = []
hc_fc.append(grid_forecast.data)
hc_fc.append(grid_hindcast.data)
hindcast_full = xr.concat(hc_fc,dim='time') ## stacking the data along month dimension
hindcast_full = hindcast_full.rename(var)

era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]
grid_observations = Observations(
                            name='Era',
                            observations=era,
                            forecast=hindcast_full,
                            process=True
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

fc_member          = list(forecast.groupby('member'))


for n,(nn,fcdata_m) in enumerate(fc_member): # loop through each member
  
  
    for lt in steps:
      
        fcdata_steps        = fcdata_m.sel(step=pd.Timedelta(lt,'D')) #loop through each month
        hc_steps            = hindcast_mean.sel(step=pd.Timedelta(lt,'D'))
        
        dim                 = 'validation_time.month'
        x_group             = list(fcdata_steps.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
        y_group             = list(hc_steps.groupby(dim))
        
        for m,(mm,xdata) in enumerate(x_group): # loop through each member
            mm_y,ydata          = y_group[m]
            
            dim                 = 'validation_time.day'
            x_group_day             = list(xdata.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
            y_group_day             = list(ydata.groupby(dim))

    
    
    fc_group_m_t          = list(fcdata_m.groupby(dim)) #loop through each month
    for nm,(nnm,fcdata_m_t) in enumerate(fc_group_m_t):
        print(nm)
        print(fcdata_m_t.dims)
  
  forecast_a        = fcdata - hindcast_mean.sel(time='2017')


## loop through each member to calculate the anomaly
  
#forecast_a        = xh.assign_validation_time(grid_forecast.data_a)
#forecast_mean     = xh.assign_validation_time(grid_forecast.mean)
#forecast_std      = xh.assign_validation_time(grid_forecast.std)



print('\tLoop through all steps')

for lt in steps:

    fc_a          = forecast_a.sel(step=pd.Timedelta(lt,'D'))
    era_a         = stacked_era_a.sel(step=pd.Timedelta(lt,'D'))
    hc_c          = hindcast_a.sel(step=pd.Timedelta(lt,'D'))
    
    
   
  
    fc_group      = list(fc_a.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    era_group     = list(era_a.groupby(dim))
    hc_group      = list(hc_c.groupby(dim))
    
   



