import pandas as pd

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
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
all_observations    = BarentsWatch().load('all',no=0).sortby('time')[var]
obs350    = BarentsWatch().load('all',no=350).sortby('time')[var]

boxes = [
            ( ( 5.5, 8.5, 57.5, 58.7 ), 'g' ),
            ( ( 4, 7.5, 58.7, 62.02 ), 'f' ),
            ( ( 4, 7.5, 6.5, 10, 62.02, 63 ), 'e' ),
            ( ( 6.5, 10, 13, 16.5, 63, 67.5), 'd' ),
            ( ( 12.5, 16.7, 67.5, 69.5), 'c' ),
            ( ( 16.7, 20.7, 16.7, 18, 67.8, 69.5, 70.3, 70.3), 'b' ),
            ( ( 20.25, 24.999, 69, 71.5), 'a' )
]

graphics.point_map(all_observations,bw_special=True,boxes=boxes,bw_poi=obs350)

exit()
















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

graphics.qq_plot(
                hisdalen_obs,
                hisdalen_ec,
                y_axlabel = 'EC',
                filename='EC'
                )

graphics.qq_plot(
                hisdalen_obs,
                hisdalen_ec.mean('member'),
                y_axlabel = 'EC',
                filename='EC_mean',
                title='ensemble mean'
                )

graphics.qq_plot(
                hisdalen_obs,
                hisdalen_pers,
                y_axlabel = 'pers',
                filename='pers'
                )

graphics.qq_plot(
                hisdalen_obs,
                hisdalen_combo,
                y_axlabel = 'combo',
                filename='combo_clustered'
                )

graphics.qq_plot(
                hisdalen_obs,
                hisdalen_combo.mean('member'),
                y_axlabel = 'combo',
                filename='combo_clustered_mean',
                title='ensemble mean'
                )
