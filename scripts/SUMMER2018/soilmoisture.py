import xarray as xr
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

from S2S.local_configuration import config





dir = '/nird/projects/NS9001K/sso102/SOILMOISTURE/'

file_monmean = 'ERA_soilmoisture_eur_2018_monmean_anom.nc'
file_ymonmean = 'ERA_soilmoisture_eur_1999_2019_ymonmean.nc'

dataopen_monmean = xr.open_dataset(dir + file_monmean)
dataopen_ymonmean = xr.open_dataset(dir + file_ymonmean)





## Plotting climatology

varplot = dataopen_ymonmean.swvl1
#levels_plot = np.linspace(round(np.nanmin(dataopen_ymonmean.swvl1)),round(np.nanmax(dataopen_ymonmean.swvl1)),21)
#levels_cbar = np.linspace(round(np.nanmin(dataopen_ymonmean.swvl1)),round(np.nanmax(dataopen_ymonmean.swvl1)),11)
levels_plot = np.linspace(0.2,0.6,21)
levels_cbar = np.linspace(0.2,0.6,11)
plot_title  = 'Climatology (1999-2019) ' + dataopen_ymonmean.swvl1.long_name 
fname       = 'Soilmoisture_swvl1'
label_text  = dataopen_ymonmean.swvl1.units

im = varplot.plot(
  x               = 'lon',
  y               = 'lat',
  col              = 'time',
  col_wrap         = 3,
  levels           = levels_plot,
  subplot_kws      = dict(projection=ccrs.PlateCarree()),
  transform        = ccrs.PlateCarree(),
  cbar_kwargs      = {'label': label_text, 'ticks': levels_cbar}
        )

  
for i,ax in enumerate(im.axes.flat):
  ax.coastlines(resolution = '10m', 
                color      = 'black',
                linewidth  = 0.2)
  ax.set_title(dataopen_ymonmean.time[i].dt.month.values)
        
plt.suptitle(plot_title)

plt.savefig(fname+'.png',dpi='figure',bbox_inches='tight')
plt.close()
print('Figure stored at: '+fname+'.png')




## Plotting 2018 anomaly

varplot = dataopen_monmean.swvl1 
levels_plot = np.linspace(round(np.nanmin(varplot)),round(np.nanmax(varplot)),21)
levels_cbar = np.linspace(round(np.nanmin(varplot)),round(np.nanmax(varplot)),11)
#levels_plot = np.linspace(0.2,0.6,21)
#levels_cbar = np.linspace(0.2,0.6,11)
plot_title  = '2018 anomaly ' + dataopen_ymonmean.swvl1.long_name
fname       = 'Soilmoisture_swvl1_2018_anom'
label_text  = dataopen_ymonmean.swvl1.units

im = varplot.plot(
  x               = 'lon',
  y               = 'lat',
  col              = 'time',
  col_wrap         = 3,
  levels           = levels_plot,
  subplot_kws      = dict(projection=ccrs.PlateCarree()),
  transform        = ccrs.PlateCarree(),
  cbar_kwargs      = {'label': label_text, 'ticks': levels_cbar}
        )

  
for i,ax in enumerate(im.axes.flat):
  ax.coastlines(resolution = '10m', 
                color      = 'black',
                linewidth  = 0.2)
  ax.set_title(dataopen_ymonmean.time[i].dt.month.values)
        
plt.suptitle(plot_title)

plt.savefig(fname+'.png',dpi='figure',bbox_inches='tight')
plt.close()
print('Figure stored at: '+fname+'.png')

