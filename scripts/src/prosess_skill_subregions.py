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
#bounds = (0,28,55,75)


high_res = False
steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

# observations must be weekly mean values with a time dimension
# load grid for catchment
data_path = '/nird/projects/NS9853K/PROJECTS/STATKRAFT/DATA/'

with open(data_path+"catchment_boundaries_wgs84.p", "rb") as stream:
    polygon_dict = pickle.load(stream)

for k in polygon_dict:
     
    minlon,minlat,maxlon,maxlat = polygon_dict[k].bounds 
    bounds = (minlon,maxlon,minlat,maxlat)

## la den lese inn fult område, også velg ut område etterpå.    
    



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
