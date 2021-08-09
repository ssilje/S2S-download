import pandas as pd
import xarray as xr

from S2S.data_handler import ERA5, ECMWF_S2SF, ECMWF_S2SH
from S2S.process      import Hindcast, Observations, Grid2Point

bounds    = (5,10,60,65) # test_domain
variables = ['u10','v10']
# variables = ['sst']

t_start  = (2020,4,23)
t_end    = (2020,4,30)

high_res = False

for var in variables:

    data = ERA5(high_res=high_res).load(
                                    var,
                                    t_start,
                                    t_end,
                                    bounds
                                )[var]
    print(data)
