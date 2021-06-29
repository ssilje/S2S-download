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

pers  = models.persistence(
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
hisdalen_combo = combo.sel(location=loc('Hisdalen'))
hisdalen_pers  = pers.sel(location=loc('Hisdalen'))
hisdalen_ec    = point_hindcast.data_a.sel(location=loc('Hisdalen'))

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            dim         = 'validation_time.month',
            filename    = 'BW_EC_vs_clim',
            title       = 'BW EC vs Clim',
            mlabs       = ['EC']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_pers],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            dim         = 'validation_time.month',
            filename    = 'BW_persistence_vs_clim',
            title       = 'BW persistence vs clim',
            mlabs       = ['pers']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = hisdalen_pers,
            dim         = 'validation_time.month',
            filename    = 'BW_EC_vs_persistence',
            title       = 'BW EC vs persistence',
            mlabs       = ['']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = hisdalen_ec,
            dim         = 'validation_time.month',
            filename    = 'BW_combo_vs_EC',
            title       = 'BW combo vs EC',
            mlabs       = ['']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = hisdalen_pers,
            dim         = 'validation_time.month',
            filename    = 'BW_combo_vs_pers',
            title       = 'BW combo vs pers',
            mlabs       = ['']
        )
