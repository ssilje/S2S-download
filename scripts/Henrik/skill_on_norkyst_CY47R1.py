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
import S2S.scoring as sc

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
bw       = BarentsWatch().load('all',no=0).sortby('time')[var]

t_start  = (2020,12,1) #can start with 8
t_end    = (2021,1,1)
model    = 'CY47R1'

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
                        process    = True,
                        download   = False,
                        split_work = False,
                        period     = [nk.time.min(),nk.time.max()]
                    )
    print([nk.time.min(),nk.time.max()])
    # update times to the next month
    t_start = t_end
    t_end   = add_month(t_end)

    observations = Observations(
                                name='NorKyst-800',
                                observations=nk,
                                forecast=hindcast,
                                process=True
                                )
    del nk

    hindcast = Grid2Point(observations,hindcast)\
                            .correlation(step_dependent=True)

    mae_fc = xs.mae(
            hindcast.data_a.mean('member'),
            observations.data_a,
            dim=[]
        )

    mae_fc = mae_fc.where(mae_fc.time.dt.month==int(month),drop=True)\
                                    .mean('time',skipna=True)

    crps_fc = sc.crps_ensemble(observations.data_a,hindcast.data_a)

    crps_fc = crps_fc.where(crps_fc.time.dt.month==int(month),drop=True)\
                                    .mean('time',skipna=True)

    crps_clim = xs.crps_gaussian(
                                    observations.data_a,
                                    xr.zeros_like(observations.data_a),
                                    xr.ones_like(observations.data_a),
                                    dim=[]
                                )

    crps_clim = crps_clim.where(crps_clim.time.dt.month==int(month),drop=True)\
                                    .mean('time',skipna=True)

    del hindcast

    pers    = models.persistence(
                    init_value   = observations.init_a,
                    observations = observations.data_a
                    )

    mae_pers = xs.mae(
            pers,
            observations.data_a,
            dim=[]
        )

    mae_pers = mae_pers.where(mae_pers.time.dt.month==int(month),drop=True)\
                                    .mean('time',skipna=True)

    del pers

    mae_clim = xs.mae(
            xr.zeros_like(observations.data_a),
            observations.data_a,
            dim=[]
        )

    mae_clim = mae_clim.where(mae_clim.time.dt.month==int(month),drop=True)\
                                    .mean('time',skipna=True)

    del observations

    for step in steps:

        for ref_fc,lab in zip([mae_clim,mae_pers],['CLIM','PERS']):

            latex.set_style(style='white')
            fig,ax = plt.subplots(1,1,\
                figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
                subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

            ss = ( 1 - mae_fc/ref_fc ).sel(step=step)

            cmap   = latex.skill_cmap().reversed()
            levels = np.arange(-1.,1.05,0.05)
            norm   = BoundaryNorm(levels,cmap.N)

            cs = ax.scatter(
                        ss.lon.values,
                        ss.lat.values,
                        c=ss.values,
                        s=1.1,
                        cmap=cmap,
                        norm=norm,
                        alpha=0.95,
                        transform=ccrs.PlateCarree()
                    )
            ax.coastlines(resolution='10m', color='grey',\
                                    linewidth=0.2)
            ax.set_title(mparser[month] + ' MAEss EC vs. '+lab+', NorKyst, lt:'\
                                                                +str(step.days))
            fig.colorbar(cs,ax=ax)
            graphics.save_fig(fig,
                            model+'mae_skill_map_NorKyst_month'+month+lab+str(step.days)
                            )

        latex.set_style(style='white')
        fig,ax = plt.subplots(1,1,\
            figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
            subplot_kw=dict(projection=ccrs.NorthPolarStereo()))

        ss = ( 1 - crps_fc/crps_clim ).sel(step=step)

        cmap   = latex.skill_cmap().reversed()
        levels = np.arange(-1.,1.05,0.05)
        norm   = BoundaryNorm(levels,cmap.N)

        cs = ax.scatter(
                    ss.lon.values,
                    ss.lat.values,
                    c=ss.values,
                    s=1.1,
                    cmap=cmap,
                    norm=norm,
                    alpha=0.95,
                    transform=ccrs.PlateCarree()
                )
        ax.coastlines(resolution='10m', color='grey',\
                                linewidth=0.2)
        ax.set_title(mparser[month] + ' CRPss EC vs. CLIM, NorKyst, lt:'\
                                                            +str(step.days))
        fig.colorbar(cs,ax=ax)
        graphics.save_fig(fig,
                        model+'crps_skill_map_NorKyst_month'+month+'CLIM'+str(step.days)
                        )
