import pandas as pd
import numpy as np
import gridpp
import json
import os

from S2S.date_helpers import get_forcast_date_cycle
from S2S.file_handling import read_grib_file
from S2S.interpolation import interpolate_and_save
from S2S.local_configuration import config
from S2S.gridpp_helpers import make_grid_object

var_name = 'sst'
var_name_abbr = 'sst'
mdl_vrsn = 'CY46R1'
S2S_dirbase = config['S2S_DIR']
product = 'forecast'

# sst
lead_time = np.arange(1, 47)
dates_fc_cycle = get_forcast_date_cycle(
    start_year=2019,
    start_month=7,
    start_day=1,
)

# sav300
# dates_fc_cycle=get_forcast_date_cycle(
#     start_year=2020,
#     start_month=5,
#     start_day=4,
# )

date_str = dates_fc_cycle[0].strftime('%Y-%m-%d')

ds_cf = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='cf',
    date_str=date_str,
)

ds_pf = read_grib_file(
    S2S_dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    cast_type='pf',
    date_str=date_str,
)

ECMWF_grid = make_grid_object(ds_cf)

with open(os.path.join(config['BW_DIR'], "metadata_BW_sites.json")) as json_file:
    data_BW = pd.DataFrame(json.load(json_file))
BW_points = gridpp.Points(data_BW.lat, data_BW.lon)

interpolate_and_save(
    var_name=var_name,
    product=product,
    date_str=date_str,
    ds_cf=ds_cf,
    ds_pf=ds_pf,
    in_grid=ECMWF_grid,
    out_points=BW_points,
    data_export_dir=config['CF_DATA_DIR'],
    data_BW=data_BW,
)

# %%
