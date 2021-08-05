import pandas as pd

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models, location_cluster

def _loc_from_name(name):
    return str(location_cluster.loc_from_name(name))

bounds   = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)


high_res = True
steps    = pd.to_timedelta([9,16,23,30,37],'D')

# observations must be weekly mean values with a time dimension
bw       = BarentsWatch().load('all',no=350).sortby('time')[var]

### get observations ###
nk = []
for loc in bw.location.values:

    print(loc)
    fname =                 config['NORKYST'] +\
                                'NorKyst800_' +\
                                     str(loc) +\
                                        '.nc'

    nk.append(xr.open_dataset(fname)[var])

nk = xr.concat(nk,'location',join='outer').drop('radius')

print(bw)
print(nk)
