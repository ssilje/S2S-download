import xarray as xr
import pandas as pd
import json
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics import latex

path   = '/nird/projects/NS9853K/DATA/norkyst800/'
fn1    = 'norkyst800_sst_'
fn2    = '_var.nc'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']

latex.set_style(style='white')
fig,axes = plt.subplots(3,4,\
    figsize=latex.set_size(width=345,subplots=(3,4),fraction=0.95),\
    subplot_kw=dict(projection=ccrs.PlateCarree()))

for ax,month in zip(axes.flatten(),months):

    print(month)

    data = xr.open_dataset(path+fn1+month+fn2).squeeze()

    print(data)
    exit()

    ax.coastlines(resolution='10m', color='grey',\
                            linewidth=0.2)

    ax.contourf(lons,lats,var,projection=ccrs.PlateCarree())

    ax.set_extent((0,25,55,75),crs=ccrs.PlateCarree())
    ax.set_title(month)
