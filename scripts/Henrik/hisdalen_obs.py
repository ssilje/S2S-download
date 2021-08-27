import xarray      as xr
import pandas      as pd
import numpy       as np
import xskillscore as xs

import json
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics            import latex, graphics

from matplotlib.colors       import BoundaryNorm

from S2S.data_handler        import ERA5, BarentsWatch
from S2S.process             import Hindcast, Observations, Grid2Point
from S2S                     import models
import S2S.scoring          as sc
import S2S.location_cluster as lc

var      = 'sst'
loc_name = 'Hisdalen'
loc      = str(lc.loc_from_name(loc_name))

fname        = config['NORKYST'] + 'NorKyst800_' + loc + '.nc'
norkyst      = xr.open_dataset(fname)[var]

barentswatch = BarentsWatch().load([loc_name]).sortby('time')[var]

norkyst = norkyst.rolling(time=7,center=True).mean()


fname = 'timeseries/timeseries_norkyst_vs_bw_'+loc_name

latex.set_style(style='white')
fig,ax = plt.subplots(1,1,
                figsize=latex.set_size(width='thesis',
                    subplots=(1,1))
                )

ax.plot(
        barentswatch.time.squeeze(),
        barentswatch.squeeze(),
        'd',
        alpha     = 0.8,
        color     = 'grey',
        ms        = 0.5,
        # linewidth = 0.5,
        label     = 'BarentsWatch'
    )

ax.plot(
        norkyst.time.squeeze(),
        norkyst.squeeze(),
        'o',
        alpha     = 0.8,
        color     = 'k',
        ms        = 0.5,
        # linewidth = 0.5,
        label     = 'Norkyst-800'
    )

ax.legend(prop={'size': 6})
fig.suptitle('Observations at '+loc_name)
ax.set_ylim(0,25)

graphics.save_fig(fig,fname)
