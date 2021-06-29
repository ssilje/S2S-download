import pandas as pd
import xarray as xr

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics
from S2S import models, location_cluster

def loc(name):
    return str(location_cluster.loc_from_name(name))

domainID = 'norwegian_coast'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)


high_res = False
steps    = pd.to_timedelta([9,16,23,30,37],'D')

# observations must be weekly mean values with a time dimension
point_observations    = BarentsWatch().load('all',no=350).sortby('time')[var]

grid_hindcast     = Hindcast(
                        var,
                        t_start,
                        t_end,
                        domainID,
                        high_res=high_res,
                        steps=steps,
                        process=False,
                        download=False,
                        split_work=True
                    )

point_observations = Observations(
                            name='BarentsWatch',
                            observations=point_observations,
                            forecast=grid_hindcast,
                            process=False
                            )

point_hindcast     = Grid2Point(point_observations,grid_hindcast)\
                            .correlation(step_dependent=True)

pers  = models.persistence(
                init_value   = point_observations.init_a,
                observations = point_observations.data_a
                )

combo = models.combo(
                        init_value      = point_observations.init_a,
                        model           = point_hindcast.data_a,
                        observations    = point_observations.data_a
                    )

combo = point_hindcast.data_a - point_hindcast.data_a.mean('member') + combo

# adjust spread
combo = models.bias_adjustment_torralba(
                            forecast        = combo,
                            observations    = point_observations.data_a,
                            spread_only     = True
                            )

mae.map(
        observations = point_observations.data_a,
        model        = point_hindcast.data_a,
        clim_mean    = xr.full_like(point_observations.data_a,0),
        dim          = 'validation_time.month',
        filename     = 'EC_vs_clim'
        )

mae.map(
        observations = point_observations.data_a,
        model        = point_hindcast.data_a,
        clim_mean    = pers,
        dim          = 'validation_time.month',
        filename     = 'EC_vs_pers'
        )

mae.map(
        observations = point_observations.data_a,
        model        = combo,
        clim_mean    = xr.full_like(point_observations.data_a,0),
        dim          = 'validation_time.month',
        filename     = 'combo_vs_clim'
        )

mae.map(
        observations = point_observations.data_a,
        model        = combo,
        clim_mean    = pers,
        dim          = 'validation_time.month',
        filename     = 'combo_vs_pers'
        )

mae.map(
        observations = point_observations.data_a,
        model        = combo,
        clim_mean    = point_hindcast.data_a,
        dim          = 'validation_time.month',
        filename     = 'combo_vs_EC'
        )
