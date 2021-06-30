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


high_res = True
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
                            process=True
                            )

point_hindcast     = Grid2Point(point_observations,grid_hindcast)\
                            .correlation(step_dependent=True)

o  = point_observations.data_a
oi = point_observations.init_a
h  = point_hindcast.data_a

obs = []
mod = []

for loc in o.location:

    o_cluster  = location_cluster.cluster(
                                        da=o,
                                        loc=str(loc.values),
                                        lon_tolerance=0.5,
                                        lat_tolerance=0.25
                                    )

    oi_cluster = oi.sel(location=o_cluster.location)
    h_cluster  = h.sel(location=o_cluster.location)



    combo = models.combo(
                            init_value      = oi_cluster,
                            model           = h_cluster,
                            observations    = o_cluster,
                            cluster_name    = 'location'
                        )

    combo = h_cluster - h_cluster.mean('member') + combo

    # adjust spread
    combo = models.bias_adjustment_torralba(
                                forecast        = combo,
                                observations    = o_cluster,
                                spread_only     = True,
                                cluster_name    = 'location'
                                )

    obs.append(o_cluster.sel(location=loc))
    mod.append(combo.sel(location=loc))

obs = xr.concat(obs,'location')
mod = xr.concat(mod,'location')

mae.map(
        observations = obs,
        model        = mod,
        clim_mean    = xr.full_like(obs,0),
        dim          = 'validation_time.month',
        filename     = 'COMBO_vs_clim'
        )

mae.map(
        observations = obs,
        model        = mod,
        clim_mean    = xr.full_like(obs,0),
        dim          = 'validation_time.season',
        filename     = 'COMBO_vs_clim'
        )
