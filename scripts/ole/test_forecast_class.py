############
# try loading hindcasts and check out the struture of the data to copy for forecasts
############

# load packages
from S2S.process import Forecast
import pandas as pd
from S2S.local_configuration import config

# define input to the hindcast class:
var      = 'sst'

t_start  = (2019,12,1)
t_end    = (2020,1,10)

high_res = False
steps    = pd.to_timedelta(range(4,46,7),'D')

bounds = (5,20,55,70)


forecast = Forecast(
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

### create some example plots to explore the data
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

figpath = config['SAVEFIG']

# plot time series (abs, anom, stdzd anoms) at one grid point and for one lead time:
# Select one location:

ltd = 18
LONIX,LATIX = 0,0

hc_loc = forecast.hc_abs.sel(step='{0:d} days'.format(ltd)).isel(lat = LATIX, lon = LONIX)
hc_mean_loc = forecast.hc_mean.sel(step='{0:d} days'.format(ltd)).isel(lat = LATIX, lon = LONIX)

# climatology:

clim_mean = forecast.hc_mean.sel(step='{0:d} days'.format(ltd),time=slice('2018-8-1','2019-7-1')).isel(lon=LONIX,lat=LATIX)
clim_std = forecast.hc_std.sel(step='{0:d} days'.format(ltd),time=slice('2018-8-1','2019-7-1')).isel(lon=LONIX,lat=LATIX)
clim_int = [clim_mean-clim_std.values,clim_mean+clim_std.values]


F,AX = plt.subplots(2,1,figsize=(7,5))
AX[0].plot(forecast.data.time.values,forecast.data.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].values.T,color='lightgrey',alpha=.5)
AX[0].plot(forecast.data.time.values,forecast.data.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].mean('member').values,color='k')

matchyear = [datetime(tt.year+1,tt.month,tt.day) for tt in clim_mean.time.values.astype('M8[ms]').astype('O')]
AX[0].plot(matchyear,clim_mean,color='r')
AX[0].fill_between(matchyear,clim_int[0],clim_int[-1],color='r',alpha=.3)

AX[0].set_title('absolute 7-day mean ssts {0:d}d lead time'.format(ltd))
AX[0].grid()
AX[0].set_xlim([forecast.data_a.time[0].values,forecast.data_a.time[-1].values])

# AX[1].plot(forecast.data_anom.time.values,forecast.data_anom.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].values.T,color='lightgrey',alpha=.5)
# AX[1].plot(forecast.data_anom.time.values,forecast.data_anom.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].mean('member').values,color='k')
# AX[1].set_title('ssts 7-day mean anomalies {0:d}d lead time'.format(ltd))
# AX[1].grid()

AX[1].plot(forecast.data_a.time.values,forecast.data_a.sel(step='{0:d} days'.format(ltd)).isel(lon=LONIX,lat=LATIX).values.T,color='lightgrey',alpha=.5)
AX[1].plot(forecast.data_a.time.values,forecast.data_a.sel(step='{0:d} days'.format(ltd)).isel(lon=LONIX,lat=LATIX).mean('member').values,color='k')

AX[1].plot(forecast.data.time.values,np.zeros_like([forecast.data.time.values],dtype=int).squeeze(),color='r')
AX[1].fill_between(forecast.data.time.values,(-1)*np.ones_like([forecast.data.time.values],dtype=int).squeeze(),
np.ones_like([forecast.data.time.values],dtype=int).squeeze(),color='r',alpha=.3)

AX[1].set_title('sst 7-day mean standardized anomalies {0:d}d lead time'.format(ltd))
AX[1].grid()
AX[1].set_xlim([forecast.data_a.time[0].values,forecast.data_a.time[-1].values])

F.subplots_adjust(hspace=.4)

F.savefig('{0:s}sst_7dm_ECMWF_fc_Jan20_lt{1:d}d_onegp.png'.format(figpath,ltd),dpi=300,bbox_inches='tight')


### try loading ERA5 matching the forecast dates

# from S2S.process import Observations
# from S2S.data_handler import ERA5

# observations = ERA5(high_res=False).load(
#                                         var=var,
#                                         start_time=t_start,
#                                         end_time=t_end,
#                                         bounds=bounds
#                                         ) - 273.15

# obs_E5 = Observations(name='ERA5',observations=observations,forecast=forecast)



# F,AX = plt.subplots(2,1,figsize=(7,5))
# AX[0].plot(forecast.data.time.values,forecast.data.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].values.T,color='lightgrey',alpha=.5)
# AX[0].plot(forecast.data.time.values,forecast.data.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].mean('member').values,color='k')

# AX[0].plot(obs_E5.data.time.values,obs_E5.data.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].values,color='C0')

# matchyear = [datetime(tt.year+1,tt.month,tt.day) for tt in clim_mean.time.values.astype('M8[ms]').astype('O')]
# AX[0].plot(matchyear,clim_mean,color='r')
# AX[0].fill_between(matchyear,clim_int[0],clim_int[-1],color='r',alpha=.3)

# AX[0].set_title('absolute 7-day mean ssts {0:d}d lead time'.format(ltd))
# AX[0].grid()
# AX[0].set_xlim([forecast.data_a.time[0].values,forecast.data_a.time[-1].values])

# # AX[1].plot(forecast.data_anom.time.values,forecast.data_anom.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].values.T,color='lightgrey',alpha=.5)
# # AX[1].plot(forecast.data_anom.time.values,forecast.data_anom.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].mean('member').values,color='k')
# # AX[1].set_title('ssts 7-day mean anomalies {0:d}d lead time'.format(ltd))
# # AX[1].grid()

# AX[1].plot(forecast.data_a.time.values,forecast.data_a.sel(step='{0:d} days'.format(ltd)).isel(lon=LONIX,lat=LATIX).values.T,color='lightgrey',alpha=.5)
# AX[1].plot(forecast.data_a.time.values,forecast.data_a.sel(step='{0:d} days'.format(ltd)).isel(lon=LONIX,lat=LATIX).mean('member').values,color='k')

# AX[1].plot(obs_E5.data_a.time.values,obs_E5.data_a.sel(step='{0:d} days'.format(ltd))[...,LONIX,LATIX].values,color='C0')

# AX[1].plot(forecast.data.time.values,np.zeros_like([forecast.data.time.values],dtype=int).squeeze(),color='r')
# AX[1].fill_between(forecast.data.time.values,(-1)*np.ones_like([forecast.data.time.values],dtype=int).squeeze(),
# np.ones_like([forecast.data.time.values],dtype=int).squeeze(),color='r',alpha=.3)

# AX[1].set_title('sst 7-day mean standardized anomalies {0:d}d lead time'.format(ltd))
# AX[1].grid()
# AX[1].set_xlim([forecast.data_a.time[0].values,forecast.data_a.time[-1].values])

# F.subplots_adjust(hspace=.4)

# F.savefig('{0:s}sst_7dm_ECMWF_ERA5_fc_Jan20_lt{1:d}d_onegp.png'.format(figpath,ltd),dpi=300,bbox_inches='tight')
