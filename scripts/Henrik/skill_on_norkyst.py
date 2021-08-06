import xarray as xr
import pandas as pd
import numpy as np
import xskillscore as xs

import json
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics import latex, graphics

from matplotlib.colors import BoundaryNorm

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point
from S2S import models

bounds = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

steps    = pd.to_timedelta([9,16,23,30,37],'D')

high_res = True

mparser = {
            '01':'JAN','02':'FEB','03':'MAR',
            '04':'APR','05':'MAY','06':'JUN',
            '07':'JUL','08':'AUG','09':'SEP',
            '10':'OCT','11':'NOV','12':'DEC'
        }
months  = ['01','02','03','04','05','06','07','08','09','10','11','12']

mae_fc, mae_clim, mae_pers = [], [], []
################################################################################
bw       = BarentsWatch().load('all',no=0).sortby('time')[var]

for loc in bw.location.values:

    print(loc)
    fname =                 config['NORKYST'] +\
                                'NorKyst800_' +\
                                     str(loc) +\
                                        '.nc'

    nk = xr.open_dataset(fname)[var].drop('radius')

    hindcast     = Hindcast(
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
    print(hindcast.data_a.step.isel(step=0).values)
    observations = Observations(
                                name='NorKyst-800_'+str(loc),
                                observations=nk,
                                forecast=hindcast,
                                process=True
                                )

    clim_fc = models.clim_fc(observations.mean,observations.std)
    pers    = models.persistence(
                    init_value   = observations.init_a,
                    observations = observations.data_a
                    )

    mae_fc.append(xs.mae(hindcast.data_a.mean('member'),observations.data_a,dim=[]))
    mae_clim.append(xs.mae(clim_fc,observations.data_a,dim=[]))
    mae_pers.append(xs.mae(pers,observations.data_a,dim=[]))

mae_fc = xr.concat(mae_fc,'location',join='outer')
mae_clim = xr.concat(mae_clim,'location',join='outer')
mae_pers = xr.concat(mae_pers,'location',join='outer')

for month in months:

    mfc   = mae_fc.where(mae_fc.time.dt.month==int(month),drop=True)\
                                                .mean('time',skipna=True)
    mclim = mae_clim.where(mae_fc.time.dt.month==int(month),drop=True)\
                                                .mean('time',skipna=True)
    mpers = mae_pers.where(mae_fc.time.dt.month==int(month),drop=True)\
                                                .mean('time',skipna=True)

    for step in steps:

        for ss,lab in zip([mfc/mclim,mfc/mpers],['CLIM','PERS']):

            latex.set_style(style='white')
            fig,ax = plt.subplots(1,1,\
                figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
                subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

            ss = ss.sel(step=step)
            ss = ( 1 - ss ).squeeze().transpose('lat','lon')

            cmap   = latex.cm_rgc(c='yellow')
            levels = np.arange(0,7,0.5)
            norm   = BoundaryNorm(levels,cmap.N)

            cs = ax.contourf(ss.lon,ss.lat,ss,transform=ccrs.PlateCarree(),
                                cmap=cmap,norm=norm,extend='max',levels=levels)
            ax.coastlines(resolution='10m', color='grey',\
                                    linewidth=0.2)
            ax.set_title(mparser[month] + ' MAEss EC vs. '+lab+', NorKyst, lt:'+str(step.values.days))
            fig.colorbar(cs,ax=ax)
            graphics.save_fig(fig,'variance_map_NorKyst_'+lab+month+str(step.values.days))
