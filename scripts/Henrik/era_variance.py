import xarray as xr
import pandas as pd
import numpy as np
import json
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics import latex, graphics

from matplotlib.colors import BoundaryNorm

from S2S.data_handler import ERA5, BarentsWatch

bounds = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

high_res = True

mparser = {
            '01':'JAN','02':'FEB','03':'MAR',
            '04':'APR','05':'MAY','06':'JUN',
            '07':'JUL','08':'AUG','09':'SEP',
            '10':'OCT','11':'NOV','12':'DEC'
        }
months  = ['01','02','03','04','05','06','07','08','09','10','11','12']

era = ERA5(high_res=True).load(
                            var         = var,
                            start_time  = clim_t_start,
                            end_time    = clim_t_end,
                            bounds      = bounds
                        )

# latex.set_style(style='white')
# fig,axes = plt.subplots(1,1,\
#     figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
#     subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

# for ax,month in zip(axes.flatten(),months):

for month in months:

    latex.set_style(style='white')
    fig,ax = plt.subplots(1,1,\
        figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
        subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

    data = era.where(era.time.dt.month==int(month),drop=True)

    data = data.sst.var('time',skipna=True).squeeze().transpose('lat','lon')

    lons = data.lon
    lats = data.lat

    cmap   = latex.cm_rgc(c='yellow')
    levels = np.arange(0,7,0.5)
    norm   = BoundaryNorm(levels,cmap.N)

    cs = ax.contourf(lons,lats,data,transform=ccrs.PlateCarree(),
                        cmap=cmap,norm=norm,extend='max',levels=levels)
    ax.coastlines(resolution='10m', color='grey',\
                            linewidth=0.2)
    ax.set_title(month)
    fig.colorbar(cs,ax=ax)
    graphics.save_fig(fig,'variance_map_era_'+month)
