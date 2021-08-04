import xarray as xr
import pandas as pd
import numpy as np
import json
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from S2S.local_configuration import config
from S2S.graphics import latex, graphics

from matplotlib.colors import BoundaryNorm

from S2S.data_handler import ERA5, BarentsWatch

bounds = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

high_res = True

mparser = {
            '01':'JAN','02':'FEB','03':'MAR',
            '04':'APR','05':'MAY','06':'JUN',
            '07':'JUL','08':'AUG','09':'SEP',
            '10':'OCT','11':'NOV','12':'DEC'
        }
months  = ['01','02','03','04','05','06','07','08','09','10','11','12']

era = ERA5(high_res=True).load(
                            var         = var,
                            start_time  = clim_t_start,
                            end_time    = clim_t_end,
                            bounds      = bounds
                        )

hindcast     = Hindcast(
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

observations = Observations(
                            name='ERA5',
                            observations=era,
                            forecast=hindcast,
                            process=True
                            )

print(hindcast)
print(observations)
