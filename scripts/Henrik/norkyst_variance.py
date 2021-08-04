import xarray as xr
import pandas as pd
import numpy as np
import json
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics import latex, graphics

from matplotlib.colors import BoundaryNorm

path   = '/nird/projects/NS9853K/DATA/norkyst800/'
# fn1    = 'norkyst800_sst_'
# fn2    = '_var.nc'
fname1 = 'norkyst800_sst_'
fname2 = '_daily_mean.nc'

months = ['01','02','03','04','05','06','07','08','09','10','11','12']

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

    fname = fname1 + '*-' + month + '-*' + fname2
    print(fname)

    try:
        # data = xr.open_dataset(path+fn1+month+fn2).squeeze()
        with xr.open_mfdataset( path + fname, parallel=True ) as data:

            print(data)
            # # load to memory
            data = data.load()
            print(data)
            data = data.temperature.var('time',skipna=True).squeeze()

            lons = data.longitude
            lats = data.latitude

            cmap   = latex.cm_rgc(c='yellow').reversed()
            levels = np.arange(0,7,0.5)
            norm   = BoundaryNorm(levels,cmap.N)

            cs = ax.contourf(lons,lats,data,transform=ccrs.PlateCarree(),
                                cmap=cmap,norm=norm)
            ax.coastlines(resolution='10m', color='grey',\
                                    linewidth=0.2)
            ax.set_title(month)
            fig.colorbar(cs,ax=ax,extend='max')
            graphics.save_fig(fig,'variance_map_norkyst_'+month)

    except (FileNotFoundError,np.core._exceptions.MemoryError) as e:
        pass
