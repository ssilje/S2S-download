import xarray as xr
import pandas as pd
import json

from S2S.local_configuration import config

from S2S.data_handler import ERA5, BarentsWatch, Archive
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics,brier
from S2S import models, location_cluster

def loc_from_name(loc):
    """
    Returns name of Barentswatch location from loc number
    """
    try:
        _ = int(loc)
        return loc

    except ValueError:
        with open(config['SITES'], 'r') as file:
            data = json.load(file)
            for line in data:
                if line['name']==loc:
                    return line["localityNo"]

def loc(name):
    return str(location_cluster.loc_from_name(name))

location_names = 'Hisdalen'

bounds   = (0,28,55,75)
var      = 'sst'

t_start  = (2020,2,1)
t_end    = (2021,2,1)

high_res = True
steps    = pd.to_timedelta([9,16,23,30,37],'D')

loc_name = 'Hisdalen'

### get observations ###
path     = '/nird/projects/NS9853K/DATA/norkyst800/'
filename = 'norkyst800_sst_*_daily_mean_at-BW.nc'

# open with dask
ds = xr.open_mfdataset( path + filename, parallel=True )
# load to memory, load() attribute
point_observations = location_cluster.cluster(
                                            da=ds.load()[var].drop('radius'),
                                            loc=loc_name,
                                            lon_tolerance=0.5,
                                            lat_tolerance=0.25
                                        )
### get hindcast ###
grid_hindcast = Hindcast(
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
                            name='NorKyst800_'+loc_name+'_cluster',
                            observations=point_observations,
                            forecast=grid_hindcast,
                            process=False
                            )

point_hindcast = Grid2Point(point_observations,grid_hindcast).correlation(
                                                        step_dependent=False
                                                            )

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
print(loc(loc_name))
his_mean       = point_observations.mean.sel(location=loc(loc_name))
his_std        = point_observations.std.sel(location=loc(loc_name))
hisdalen_obs   = point_observations.data_a.sel(location=loc(loc_name))
hisdalen_pers  = pers.sel(location=loc(loc_name))
hisdalen_combo = combo.sel(location=loc(loc_name))
hisdalen_ec    = point_hindcast.data_a.sel(location=loc(loc_name))

graphics.timeseries(
                        observations    = hisdalen_obs * his_std + his_mean,
                        cast            = [
                                            hisdalen_pers * his_std + his_mean,
                                            hisdalen_ec * his_std + his_mean,
                                            hisdalen_combo * his_std + his_mean
                                        ],
                        lead_time       = [9],
                        clabs           = ['PERS','EC','COMBO'],
                        filename        = 'NorKyst_persistence_combo_clustered',
                        title           = 'Obs: NorKyst800'
                    )

graphics.timeseries(
                        observations    = hisdalen_obs * his_std + his_mean,
                        cast            = [
                                            hisdalen_pers * his_std + his_mean,
                                            hisdalen_ec * his_std + his_mean,
                                            hisdalen_combo * his_std + his_mean
                                        ],
                        lead_time       = [16],
                        clabs           = ['PERS','EC','COMBO'],
                        filename        = 'NorKyst_persistence_combo_clustered',
                        title           = 'Obs: NorKyst800'
                    )

graphics.timeseries(
                        observations    = hisdalen_obs * his_std + his_mean,
                        cast            = [
                                            hisdalen_pers * his_std + his_mean,
                                            hisdalen_ec * his_std + his_mean,
                                            hisdalen_combo * his_std + his_mean
                                        ],
                        lead_time       = [23],
                        clabs           = ['PERS','EC','COMBO'],
                        filename        = 'NorKyst_persistence_combo_clustered',
                        title           = 'Obs: NorKyst800'
                    )

graphics.timeseries(
                        observations    = hisdalen_obs * his_std + his_mean,
                        cast            = [
                                            hisdalen_pers * his_std + his_mean,
                                            hisdalen_ec * his_std + his_mean,
                                            hisdalen_combo * his_std + his_mean
                                        ],
                        lead_time       = [30],
                        clabs           = ['PERS','EC','COMBO'],
                        filename        = 'NorKyst_persistence_combo_clustered',
                        title           = 'Obs: NorKyst800'
                    )

graphics.timeseries(
                        observations    = hisdalen_obs * his_std + his_mean,
                        cast            = [
                                            hisdalen_pers * his_std + his_mean,
                                            hisdalen_ec * his_std + his_mean,
                                            hisdalen_combo * his_std + his_mean
                                        ],
                        lead_time       = [37],
                        clabs           = ['PERS','EC','COMBO'],
                        filename        = 'NorKyst_persistence_combo_clustered',
                        title           = 'Obs: NorKyst800'
                    )
exit()
###############################################################################
mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            dim         = 'time.month',
            filename    = 'Norkyst_EC_vs_clim',
            title       = 'Norkyst EC vs Clim',
            mlabs       = ['EC']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_pers],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            dim         = 'time.month',
            filename    = 'Norkyst_persistence_vs_clim',
            title       = 'Norkyst persistence vs clim',
            mlabs       = ['pers']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = hisdalen_pers,
            dim         = 'time.month',
            filename    = 'Norkyst_EC_vs_persistence',
            title       = 'Norkyst EC vs persistence',
            mlabs       = ['']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = hisdalen_ec,
            dim         = 'time.month',
            filename    = 'Norkyst_combo_vs_EC',
            title       = 'Norkyst combo vs EC',
            mlabs       = ['']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = hisdalen_pers,
            dim         = 'time.month',
            filename    = 'Norkyst_combo_vs_pers',
            title       = 'Norkyst combo vs pers',
            mlabs       = ['']
        )



mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            dim         = 'time.season',
            filename    = 'Norkyst_EC_vs_clim',
            title       = 'Norkyst EC vs Clim',
            mlabs       = ['EC']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_pers],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            dim         = 'time.season',
            filename    = 'Norkyst_persistence_vs_clim',
            title       = 'Norkyst persistence vs clim',
            mlabs       = ['pers']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = hisdalen_pers,
            dim         = 'time.season',
            filename    = 'Norkyst_EC_vs_persistence',
            title       = 'Norkyst EC vs persistence',
            mlabs       = ['']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = hisdalen_ec,
            dim         = 'time.season',
            filename    = 'Norkyst_combo_vs_EC',
            title       = 'Norkyst combo vs EC',
            mlabs       = ['']
        )

mae.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = hisdalen_pers,
            dim         = 'time.season',
            filename    = 'Norkyst_combo_vs_pers',
            title       = 'Norkyst combo vs pers',
            mlabs       = ['']
        )

exit()
###########################################
std = 1.

brier.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            clim_std    = xr.full_like(hisdalen_obs,1),
            std         = std,
            dim         = 'validation_time.month',
            filename    = 'Norkyst_EC_vs_clim',
            title       = 'Norkyst EC vs Clim',
            mlabs       = ['EC']
        )


brier.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_ec],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            clim_std    = xr.full_like(hisdalen_obs,1),
            std         = std,
            dim         = 'validation_time.season',
            filename    = 'Norkyst_EC_vs_clim',
            title       = 'Norkyst EC vs Clim',
            mlabs       = ['EC']
        )

brier.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            clim_std    = xr.full_like(hisdalen_obs,1),
            std         = std,
            dim         = 'validation_time.month',
            filename    = 'Norkyst_combo_vs_clim',
            title       = 'Norkyst combo vs clim',
            mlabs       = ['']
        )

brier.skill_agg(
            in_obs      = hisdalen_obs,
            in_mod      = [hisdalen_combo],
            clim_mean   = xr.full_like(hisdalen_obs,0),
            clim_std    = xr.full_like(hisdalen_obs,1),
            std         = std,
            dim         = 'validation_time.season',
            filename    = 'Norkyst_combo_vs_clim',
            title       = 'Norkyst combo vs clim',
            mlabs       = ['']
        )
