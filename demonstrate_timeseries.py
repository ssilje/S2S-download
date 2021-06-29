import pandas as pd

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models

domainID = 'norwegian_coast'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)


high_res = True
steps    = pd.to_timedelta([9,16,23,30,37],'D')

grid_hindcast = Hindcast(
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

# observations must be weekly mean values with a time dimension
point_observations = BarentsWatch().load(
                                    [
                                        'Hisdalen',
                                        'Stokkvika'
                                    ]
                                )[var]

point_observations = Observations(
                            name='BarentsWatch',
                            observations=point_observations,
                            forecast=grid_hindcast,
                            process=True
                            )

point_hindcast = Grid2Point(point_observations,grid_hindcast).correlation(
                                                        step_dependent=True
                                                            )

clim_fc = models.clim_fc(point_observations.mean,point_observations.std)

# In absolute values
graphics.timeseries(
                        observations    = point_observations.data,
                        cast            = [clim_fc,point_hindcast.data],
                        lead_time       = [9,16],
                        clabs           = ['clim','EC'],
                        filename        = 'BW_absolute',
                        title           = 'Barentswatch EC'
                    )

# As anomalies
graphics.timeseries(
                        observations    = point_observations.data_a,
                        cast            = [point_hindcast.data_a],
                        lead_time       = [9,16],
                        clabs           = ['EC'],
                        filename        = 'BW_anomalies',
                        title           = 'Barentswatch EC'
                    )

# Simple bias adjustment
simple_bias_adjustment = ( point_hindcast.data_a * point_observations.std )\
                                + point_observations.mean
graphics.timeseries(
                        observations    = point_observations.data,
                        cast            = [clim_fc,simple_bias_adjustment],
                        lead_time       = [9,16],
                        clabs           = ['EC'],
                        filename        = 'BW_simple_bias_adjustment',
                        title           = 'Barentswatch EC'
                    )

# A slightly more complicated bias adjustment
torralba360 = models.bias_adjustment_torralba(
            forecast=point_hindcast.data-point_hindcast.mean,
            observations=point_observations.data-point_observations.mean,
            window=360
        ) + point_observations.mean

torralba30 = models.bias_adjustment_torralba(
            forecast=point_hindcast.data-point_hindcast.mean,
            observations=point_observations.data-point_observations.mean,
            window=30
        ) + point_observations.mean

graphics.timeseries(
                        observations    = point_observations.data,
                        cast            = [
                                            simple_bias_adjustment,
                                            torralba30,
                                            torralba360
                                        ],
                        lead_time       = [9,16],
                        clabs           = ['EC sb','Torralba30','Torralba360'],
                        filename        = 'BW_torralba',
                        title           = 'Barentswatch Torralba'
                    )
