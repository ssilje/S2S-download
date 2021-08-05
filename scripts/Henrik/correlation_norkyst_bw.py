import pandas as pd
import xarray as xr
import numpy as np
import json

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics,latex
from S2S import models, location_cluster

from S2S.local_configuration import config

def corr(x,y):
    """
    """
    idx_bool = np.logical_and(
                        np.isfinite(x),
                        np.isfinite(y)
                    )
    if idx_bool.sum()<2:
        r = np.nan

    else:
        r,p = stats.pearsonr(
                            x[idx_bool],
                            y[idx_bool]
                        )
    return r

rho.append(np.array(yrho))

    return np.stack(rho,axis=-1)

def _loc_from_name(name):
    return str(location_cluster.loc_from_name(name))

bounds   = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)


high_res = True
steps    = pd.to_timedelta([9,16,23,30,37],'D')

# observations must be weekly mean values with a time dimension
bw       = BarentsWatch().load('all',no=350).sortby('time')[var]

### get observations ###
nk = []
for loc in bw.location.values:

    print(loc)
    fname =                 config['NORKYST'] +\
                                'NorKyst800_' +\
                                     str(loc) +\
                                        '.nc'

    nk.append(xr.open_dataset(fname)[var])

nk = xr.concat(nk,'location',join='outer').drop('radius')

mparser = {
            '01':'JAN','02':'FEB','03':'MAR',
            '04':'APR','05':'MAY','06':'JUN',
            '07':'JUL','08':'AUG','09':'SEP',
            '10':'OCT','11':'NOV','12':'DEC'
        }
months  = ['01','02','03','04','05','06','07','08','09','10','11','12']

for month in months:

    latex.set_style(style='white')
    fig,ax = plt.subplots(1,1,\
        figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
        subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

    mbw = bw.where(bw.time.dt.month==int(month),drop=True)
    mnk = nk.where(nk.time.dt.month==int(month),drop=True)

    rho = xr.apply_ufunc(
                corr,mbw,mnk,
                input_core_dims  = [['time'],['time']],
                output_core_dims = [[]],
                vectorize=True,dask='parallelized'
                )

    rho = rho.squeeze().transpose('lat','lon')

    cmap   = latex.cm_rgc(c='white').reversed()
    levels = np.arange(-1,1.1,0.1)
    norm   = BoundaryNorm(levels,cmap.N)

    cs = ax.contourf(rho.lon,rho.lat,rho,transform=ccrs.PlateCarree(),
                        cmap=cmap,norm=norm,extend='max',levels=levels)
    ax.coastlines(resolution='10m', color='grey',\
                            linewidth=0.2)
    ax.set_title(mparser[month] + 'corr(Norkyst,Barentswatch)')
    fig.colorbar(cs,ax=ax)
    graphics.save_fig(fig,'corr_NK_BW_'+month)
