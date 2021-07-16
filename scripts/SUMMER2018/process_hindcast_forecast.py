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

steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')
#steps           = pd.to_timedelta(np.linspace(1,46,46),'D')

print('\tProcess Forecast')
#grid_forecast = Forecast(
#                        var,
#                        t_start,
#                        t_end,
#                        bounds,
#                        high_res=high_res,
#                        steps=steps,
#                        process=False,
#                        download=False,
#                        split_work=True
#                    )

grid_forecast = Forecast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        process=True,
                        download=False,
                        split_work=True
                    )



print('\tProcess Hindcast')
#grid_hindcast = Hindcast(
#                        var,
#                        t_start,
#                        t_end,
#                        bounds,
#                        high_res=high_res,
#                        steps=steps,
#                        process=False,
#                        download=False,
#                        split_work=True
#                    )

grid_hindcast = Hindcast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        process=True,
                        download=False,
                        split_work=True
                    )

## Sjekk - dataene blir berre nan. loopen under ser ok ut. 
hc_fc = []
hc_fc.append(grid_hindcast.data)
hc_fc.append(grid_forecast.data)

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



reanalysis         = xh.assign_validation_time(grid_observations.data)


hindcast           = xh.assign_validation_time(grid_hindcast.data)


forecast           = xh.assign_validation_time(grid_forecast.data)

#fc_member          = list(forecast.groupby('member'))


#for n,(nn,fcdata_m) in enumerate(fc_member): # loop through each member
fcc_step = []    
for lt in steps:
      fc_steps          = forecast.sel(step=pd.Timedelta(lt,'D')) #loop through each month
      hc_steps          = hindcast.sel(step=pd.Timedelta(lt,'D'))
        
      dim               = 'validation_time.month'
      fc_group          = list(fc_steps.groupby(dim)) 
      hc_group          = list(hc_steps.groupby(dim))
      fcc_month = []  
      for m,(mf,fcdata) in enumerate(fc_group): #loop through each month
          mh,hcdata          = hc_group[m]
            
          dim                 = 'validation_time.day'
          fc_group_day        = list(fcdata.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
          hc_group_day        = list(hcdata.groupby(dim))
          
          fcc_day = []
          for mm,(mmf,fcdata_m) in enumerate(fc_group_day): #loop through each month
              mmh,hcdata_m          = hc_group_day[mm]
             
              fcc = []
              for ensm in range(0,51,1):
                  fc_anom = hcdata_m.mean('time').mean('member') - fcdata_m.sel(member=ensm)
                  
                  fcc.append(fc_anom)
              fcc_member = xr.concat(fcc,dim='member')    
              fcc_day.append(fcc_member)
          fcc_member_day = xr.concat(fcc_day,dim='time')
          fcc_month.append(fcc_member_day)
      
      fcc_member_day_month = xr.concat(fcc_month,dim='time')
      fcc_step.append(fcc_member_day_month)
  forecast_anom = xr.concat(ffcc_step,dim='step')    
      
                 # fc_anom=fc_anom.assign_coords(time_month=xlabel)
          

    
    
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
    
   



