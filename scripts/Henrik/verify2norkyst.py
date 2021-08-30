import xarray as xr
import pandas as pd
import json

from S2S.local_configuration import config

from S2S.data_handler import ERA5, BarentsWatch, Archive
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models

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

bounds   = (0,28,55,75)
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

high_res = True
steps    = pd.to_timedelta([9,16,23,30,37],'D')

### get observations ###
path     = '/nird/projects/NS9853K/DATA/norkyst800/'
filename = 'norkyst800_sst_*_daily_mean_at-BW.nc'

location_names = ['Hisdalen','Stokkvika']

data = []
for loc_name in location_names:

    print(loc_name)
    fname =                 config['NORKYST'] +\
                                'NorKyst800_' +\
                 str(loc_from_name(loc_name)) +\
                                        '.nc'

    data.append(xr.open_dataset(fname)[var])

point_observations = xr.concat(data,'location',join='outer').drop('radius')

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
                            name='NorKyst800',
                            observations=point_observations,
                            forecast=grid_hindcast,
                            process=True
                            )

point_hindcast = Grid2Point(point_observations,grid_hindcast).correlation(
            ['Hisdalen','Stokkvika']                                            step_dependent=True
                                                            )

clim_fc = models.clim_fc(point_observations.mean,point_observations.std)

# # In absolute values
# graphics.timeseries(
#                         observations    = point_observations.data,
#                         cast            = [clim_fc,point_hindcast.data],
#                         lead_time       = [9,16],
#                         clabs           = ['clim','EC'],
#                         filename        = 'NorKystAbs_absolute',
#                         title           = 'NorKyst800 EC'
#                     )
#
# # As anomalies
# graphics.timeseries(
#                         observations    = point_observations.data_a,
#                         cast            = [point_hindcast.data_a],
#                         lead_time       = [9,16],
#                         clabs           = ['EC'],
#                         filename        = 'NorKystAnom',
#                         title           = 'NorKyst800 EC'
#                     )

# Simple bias adjustment
simple_bias_adjustment = ( point_hindcast.data_a * point_observations.std )\
                                + point_observations.mean
# graphics.timeseries(
#                         observations    = point_observations.data,
#                         cast            = [clim_fc,simple_bias_adjustment],
#                         lead_time       = [9,16],
#                         clabs           = ['NorKyst800','EC'],
#                         filename        = 'NorKyst_simple_bias_adjustment',
#                         title           = 'NorKyst800 EC'
#                     )
#
# graphics.timeseries(
#                         observations    = point_observations.data,
#                         cast            = [clim_fc,simple_bias_adjustment],
#                         lead_time       = [23,30,37],
#                         clabs           = ['NorKyst800','EC'],
#                         filename        = 'NorKyst_simple_bias_adjustment',
#                         title           = 'NorKyst800 EC'
#                     )

# more advanced modeling
clim_fc = models.clim_fc(point_observations.mean,point_observations.std)
pers    = models.persistence(
                init_value   = point_observations.init_a,
                observations = point_observations.data_a
                )

combo = models.combo(
                        init_value      = point_observations.init_a,
                        model           = point_hindcast.data_a,
                        observations    = point_observations.data_a
                    )

combo = point_hindcast.data_a - point_hindcast.data_a.mean('member') + combo

# adjust spread
combo = models.bias_adjustment_torralba(
                            forecast        = combo,
                            observations    = point_observations.data_a,
                            spread_only     = True
                            )
graphics.timeseries(
                        observations    = point_observations.data_a,
                        cast            = [pers,point_hindcast.data_a,combo],
                        lead_time       = [9,16],
                        clabs           = ['PERS','EC','combo'],
                        filename        = 'NorKyst_persistence_combo',
                        title           = 'NorKyst800 EC'
                    )

graphics.timeseries(
                        observations    = point_observations.data_a,
                        cast            = [pers,point_hindcast.data_a,combo],
                        lead_time       = [23,30,37],
                        clabs           = ['PERS','EC','combo'],
                        filename        = 'NorKyst_persistence_combo',
                        title           = 'NorKyst800 EC'
                    )
##################################################################
# start_date = '2012-06-27'
# end_date   = '2019-02-26'
#
# date_range = pd.date_range(start=start_date,end=end_date,freq='D')
#
# for n,date in enumerate(date_range):
#
#     filename_nird = 'norkyst800_sst_'+date.strftime('%Y-%m-%d')+\
#                         '_daily_mean_at-BW-loc.nc'
#
#     filename_out = 'norkyst800_sst_'+date.strftime('%Y-%m-%d')+\
#                         '_daily_mean_at-BW.nc'
#
#     print(filename_nird)
#
#     try:
#
#         ds = xr.open_dataset( path + filename_nird )
#         ds = ds.assign_coords(time=[pd.Timestamp(date)])
#         ds = ds.to_array().rename('sst').isel(variable=0).drop('variable')
#
#         ds.to_netcdf(path+filename_out)
#         print('\n')
#
#     except FileNotFoundError:
#         print('Not found\n')
#         pass
#
# def name_from_loc(loc):
#     with open(config['SITES'], 'r') as file:
#         data = json.load(file)
#         for line in data:
#             if line["localityNo"]==int(loc):
#                 return line['name']
# # open with dask
# ds = xr.open_mfdataset( path + filename, parallel=True )
# # load to memory
# da = ds.load()[var]
# # make dir
# Archive().make_dir(config['NORKYST'])
#
# for loc in da.location:
#
#     fname = 'NorKyst800_'+str(loc.values)+'.nc'
#
#     file = da.sel(location=loc).sortby('time')
#
#     file.to_netcdf(config['NORKYST'] + fname)
#
#     print(name_from_loc(loc.values))
#     print(config['NORKYST'] + fname)
