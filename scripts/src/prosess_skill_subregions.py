import pandas as pd
import xarray as xr
import pickle

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics
from S2S import models, location_cluster

def loc(name):
    return str(location_cluster.loc_from_name(name))

var             = 't2m'
clabel          = 'K'
#var             = 'abs_wind'
t_start         = (2019,7,1)
t_end           = (2020,6,26)


clim_t_start    = (1999,1,1)
clim_t_end      = (2021,1,4)
bounds = (0,28,55,75)


high_res = False
steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

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

era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]

grid_observations = Observations(
                            name='Era',
                            observations=era,
                            forecast=grid_hindcast,
                            process=False
                            )


## woring with anomalies
hindcast    = grid_hindcast.data_a
obs         = grid_observations.data_a

# load grid for catchment
data_path = '/nird/projects/NS9853K/PROJECTS/STATKRAFT/DATA/'

with open(data_path+"catchment_boundaries_wgs84.p", "rb") as stream:
    polygon_dict = pickle.load(stream)

for k in polygon_dict:
    print(k) 
    
    minlon,minlat,maxlon,maxlat = polygon_dict[k].bounds 
    bounds                      = (round(minlon-0.5),round(maxlon+0.5),round(minlat-0.5),round(maxlat+0.5))
    
    hindcast_sel    = hindcast.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
    hindcast_sel    = hindcast_sel.assign_coords(location=k)
    hindcast_sel    = hindcast_sel.mean('lon').mean('lat')
    
    obs_sel         = obs.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
    obs_sel         = obs_sel.assign_coords(location=k)
    obs_sel         = obs_sel.mean('lon').mean('lat')
    
    
    
