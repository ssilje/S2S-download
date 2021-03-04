import xarray as xr
import os
import gridpp
import numpy as np

def read_grib_file(
    S2S_dirbase,
    product,
    model_version,
    var_name_abbr,
    cast_type,
    date_str,
):
    file_name =  '_'.join([var_name_abbr, model_version, date_str, cast_type, product]) + '.grb'
    file_path = os.path.join(S2S_dirbase, product, 'ECMWF', 'sfc', var_name_abbr, file_name)
    
    print('reading file: ' + file_path)    
    dataopen = xr.open_dataset(file_path, engine='cfgrib')

    return dataopen

def make_grid(lats, lons):
    latlats, lonlons = np.meshgrid(
        lats, lons
    )
    grid = gridpp.Grid(latlats, lonlons)
    return grid