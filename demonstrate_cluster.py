import pandas as pd

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models, location_cluster

def loc(name):
    return str(location_cluster.loc_from_name(name))

bounds = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)


high_res = True
steps    = pd.to_timedelta([9,16,23,30,37],'D')

# observations must be weekly mean values with a time dimension
all_observations    = BarentsWatch().load('all',no=350).sortby('time')[var]

point_observations  = location_cluster.cluster(
                                                da=all_observations,
                                                loc='Hisdalen',
                                                lon_tolerance=0.5,
                                                lat_tolerance=0.25
                                            )

grid_hindcast     = Hindcast(
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

point_observations = Observations(
                            name='BarentsWatch',
                            observations=point_observations,
                            forecast=grid_hindcast,
                            process=True
                            )

point_hindcast     = Grid2Point(point_observations,grid_hindcast)\
                            .correlation(step_dependent=True)

clim_fc = models.clim_fc(point_observations.mean,point_observations.std)
pers    = models.persistence(
                init_value   = point_observations.init_a,
                observations = point_observations.data_a
                )

combo = models.combo(
                        init_value      = point_observations.init_a,
                        model           = point_hindcast.data_a,
                        observations    = point_observations.data_a,
                        cluster_name    = 'location'
                    )

combo = point_hindcast.data_a - point_hindcast.data_a.mean('member') + combo

# adjust spread
combo = models.bias_adjustment_torralba(
                            forecast        = combo,
                            observations    = point_observations.data_a,
                            spread_only     = True,
                            cluster_name    = 'location'
                            )

hisdalen_obs   = point_observations.data_a.sel(location=loc('Hisdalen'))
hisdalen_pers  = pers.sel(location=loc('Hisdalen'))
hisdalen_combo = combo.sel(location=loc('Hisdalen'))
hisdalen_ec    = point_hindcast.data_a.sel(location=loc('Hisdalen'))

graphics.timeseries(
                        observations    = hisdalen_obs,
                        cast            = [
                                            hisdalen_pers,
                                            hisdalen_ec,
                                            hisdalen_combo
                                        ],
                        lead_time       = [9,16],
                        clabs           = ['persistence','EC','combo'],
                        filename        = 'BW_persistence_combo_clustered',
                        title           = 'Barentswatch'
                    )

graphics.timeseries(
                        observations    = hisdalen_obs,
                        cast            = [
                                            hisdalen_pers,
                                            hisdalen_ec,
                                            hisdalen_combo
                                        ],
                        lead_time       = [23,30],
                        clabs           = ['persistence','EC','combo'],
                        filename        = 'BW_persistence_combo_clustered',
                        title           = 'Barentswatch'
                    )

graphics.timeseries(
                        observations    = hisdalen_obs,
                        cast            = [
                                            hisdalen_pers,
                                            hisdalen_ec,
                                            hisdalen_combo
                                        ],
                        lead_time       = [37],
                        clabs           = ['persistence','EC','combo'],
                        filename        = 'BW_persistence_combo_clustered',
                        title           = 'Barentswatch'
                    )
