import xarray as xr
import pandas as pd
import numpy as np
import xskillscore as xs

import json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics import latex, graphics

from matplotlib.colors import BoundaryNorm

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point
from S2S import models
import S2S.scoring as sc
import S2S.location_cluster as lc

def plus_minus_15_days(t_start,t_end):
    if t_start[1]==1:
        t_start = (t_start[0]-1,12,15)
    else:
        t_start = (t_start[0],t_start[1]-1,15)

    if t_end[1]==12:
        t_end = (t_end[0]+1,1,15)
    else:
        t_end = (t_end[0],t_end[1]+1,15)

    return t_start,t_end

# bounds = (0,28,55,75)
bounds = (0,28,55,75)
var      = 'sst'

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

steps    = pd.to_timedelta([9,16,23,30,37],'D')

high_res = True

mparser = {
            '1':'JAN','2':'FEB','3':'MAR',
            '4':'APR','5':'MAY','6':'JUN',
            '7':'JUL','8':'AUG','9':'SEP',
            '10':'OCT','11':'NOV','12':'DEC'
        }
months  = ['01','02','03','04','05','06','07','08','09','10','11','12']

mae_fc, mae_clim, mae_pers = [], [], []
################################################################################
bw = BarentsWatch().load('all',no=0).sortby('time')[var]
bw = lc.cluster(bw,'Hisdalen',3,1)

t_start  = (2020,8,1) #can start with 8
t_end    = (2021,9,1)
model    = 'CY47R1'
extent   = [4.5,7.5,59.3,61]

hh = Hindcast(
                        var,
                        (2020,1,23),
                        (2020,1,24),
                        bounds,
                        high_res=high_res,
                        steps=steps,
                        process=False,
                        download=False,
                        split_work=False,
                    )
add_month    = hh.add_month
smaller_than = hh.smaller_than
del hh

# assign new end time, one month after start time
t_end = add_month(t_start)

# until end time is reached; load one month at the time
while smaller_than(t_end,(2021,4,1)):

    month = str(t_start[1])
    print(month)

    nk = []
    for loc in bw.location.values:

        fname =                 config['NORKYST'] +\
                                    'NorKyst800_' +\
                                         str(loc) +\
                                            '.nc'

        # print(graphics.name_from_loc(loc),loc,type(loc))

        nk.append(xr.open_dataset(fname)[var].drop('radius'))
    nk = xr.concat(nk,'location')

    ts,te    = plus_minus_15_days(t_start,t_end)
    hindcast = Hindcast(
                        var,
                        ts,
                        te,
                        bounds,
                        high_res   = high_res,
                        steps      = steps,
                        process    = False,
                        download   = False,
                        split_work = False,
                        period     = [nk.time.min(),nk.time.max()]
                    )

    # update times to the next month
    t_start = t_end
    t_end   = add_month(t_end)

    observations = Observations(
                                name='Hardanger_NorKyst-800_'+str(month),
                                observations=nk,
                                forecast=hindcast,
                                process=False
                                )
    del nk
    del hindcast

    ref_loc = str(lc.loc_from_name('Klungsholmen'))
    ref_obs = observations.data_a.sel(location=ref_loc,step=steps[0]).squeeze()

    ref_obs = ref_obs.where(ref_obs.time.dt.month==int(month))
    all_obs = observations.data_a.where(
                                observations.data_a.time.dt.month==int(month)
                                )
    print(ref_obs)
    poi = ref_obs

    ref_obs = ref_obs.broadcast_like(all_obs)

    all_rho = xr.corr(ref_obs,all_obs,dim=['time'])

    all_rho = all_rho.assign_coords(location=all_rho.location.values.astype(int))
    all_obs = all_obs.assign_coords(location=all_obs.location.values.astype(int))

    all_rho = all_rho.sortby('location')
    all_obs = all_obs.sortby('location')

    all_rho = all_rho.assign_coords(location=all_obs.location)

    del observations

    for step in steps:

        latex.set_style(style='white')
        fig,ax = plt.subplots(1,1,\
            figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
            subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

        rho = all_rho.sel(step=step)

        cmap   = latex.skill_cmap()
        levels = np.arange(-1.,1.1,0.1)
        norm   = BoundaryNorm(levels,cmap.N)

        cs = ax.scatter(
                    rho.lon.values,
                    rho.lat.values,
                    c=rho.values,
                    s=3.1,
                    cmap=cmap,
                    norm=norm,
                    alpha=0.95,
                    transform=ccrs.PlateCarree(),
                    zorder=30,
                    edgecolors='k',
                    linewidth=0.2
                )

        ax.scatter(
                    poi.lon.values,
                    poi.lat.values,
                    c='k',
                    s=5,
                    alpha=1,
                    transform=ccrs.PlateCarree(),
                    zorder=40,
                    edgecolors='k',
                    linewidth=0.2
                )

        ax.set_extent(extent, crs=ccrs.PlateCarree())

        resol = '10m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
            scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land'])
        # ocean = cfeature.NaturalEarthFeature('physical', 'ocean', \
        #     scale=resol, edgecolor='none', facecolor=cfeature.COLORS['water'])
        ax.add_feature(land, zorder=1, facecolor='beige', linewidth=0.2)
        # ax.add_feature(ocean, linewidth=0.2, zorder=0 )

        # ax.coastlines(resolution='10m', color='grey',\
        #                         linewidth=0.2)
        ax.set_title(mparser[month] + 'Correlation lag: '\
                                +str(step.days-4-5)+'-'+str(step.days+3-12))
        cbar=fig.colorbar(cs,ax=ax)
        cbar.set_ticks(levels)
        cbar.set_ticklabels(levels)
        graphics.save_fig(fig,
                        model+'hardanger_correlation_NorKyst_month'+month+'_'+str(step.days)
                        )
