import xarray as xr
import numpy as np
import pandas as pd
import cartopy as cr
import matplotlib.pyplot as plt
import json

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

# local dependency
from .nird2local import pull, read
from .grid import cap_grid

def mark_box(ax,bounds,proj,col='gray',lwd = 0.75):
    """
    Mark subplot with box around bounds
    """
    lines = bounds
    ax.plot([lines[0],lines[0]],[lines[2],lines[3]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)
    ax.plot([lines[1],lines[1]],[lines[2],lines[3]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)
    ax.plot([lines[0],lines[1]],[lines[2],lines[2]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)
    ax.plot([lines[0],lines[1]],[lines[3],lines[3]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)

# load ERA-5 data
filename = '/sea_surface_temperature_era_2010_7.nc'
uni_path = '/SFE/ERA_monthly_nc'
username = 'heau'

data = read(uni_path,filename,username)

lon  = data['longitude']
lat  = data['latitude']
time = data['time']
sst  = data['sst']

# map bounds
bounds = (2,10,55,64)
lon,lat,sst = cap_grid(lon,lat,sst,bounds)

sst[sst==-32767.] = np.nan

# domain bounds
with open('./data/EIDE/domains.json') as file:
    domains = json.load(file)
domain_bounds = domains['NVS']['bounds']

# plot
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(color='black', linewidth=0.2)
cs = ax.contourf(lon,lat,sst.mean(axis=0),projection=ccrs.PlateCarree())
mark_box(ax,domain_bounds,ccrs.PlateCarree(),lwd = 1.)
plt.colorbar(cs)
plt.show()
plt.close()
