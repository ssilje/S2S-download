import pandas as pd
import numpy as np
import gridpp
import json

from S2S.date_helpers import get_forcast_date_cycle
from S2S.gridpp_helpers import make_grid_object
from S2S.file_handling import read_grib_file
from local_configuration import config

dates_fcycle = get_forcast_date_cycle(
    start_year=2020,
    start_month=5,
    start_day=4,
    num_periods=2,
)

var_name_abbr = 'sst'
mdl_vrsn = 'CY46R1'
S2S_dirbase = config['S2S_DIR']
product = 'forecast'
curr_date = dates_fcycle[0].strftime('%Y-%m-%d')

def make_points_flatten(lats, lons, valid_points):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_points = gridpp.Points(latlats.flatten()[valid_points], lonlons.flatten()[valid_points])
    return ECMWF_points


def make_grid_flatten(lats, lons, valid_points):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    ECMWF_points = gridpp.Grid(latlats.flatten()[valid_points], lonlons.flatten()[valid_points])
    return ECMWF_points

SST_GRID1_5deg = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=curr_date,
)

with open("metadata_BW_sites.json") as json_file:
    data_BW = pd.DataFrame(json.load(json_file))

BW_grid = gridpp.Points(data_BW.lat, data_BW.lon)

ECMWF_grid1_5deg = make_grid_object(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data)

SST_1_5deg2point = gridpp.nearest(
    ECMWF_grid1_5deg,
    BW_grid,
    gridpp.fill_missing(np.transpose(SST_GRID1_5deg.sst[0, :, :].data))
)
print('test 1 med grid2points')
print(SST_1_5deg2point)

sst = SST_GRID1_5deg.sst[0, :, :].data.flatten()
valid_points = np.isnan(sst) == 0  # Ocean points

ECMWF_points = make_points_flatten(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data, valid_points)

SST_points = gridpp.nearest(ECMWF_points, BW_grid, sst[valid_points])
print('test 2 med points2points')
print(SST_points)
#
# ECMWF_grid1_5deg_valid_points = make_grid_flatten(SST_GRID1_5deg.latitude.data, SST_GRID1_5deg.longitude.data,
#                                                   valid_points)
# SST_grid_points = gridpp.nearest(ECMWF_grid1_5deg_valid_points, BW_grid, sst[valid_points]