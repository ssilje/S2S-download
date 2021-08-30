import numpy as np
import xarray as xr
import pandas as pd
import json
from S2S.local_configuration import config

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

def cluster(da,loc,lon_tolerance,lat_tolerance):

    if not isinstance(loc,int):
        loc = loc_from_name(loc)

    z = da.sel(location=str(loc))

    lon_lim = np.array([-lon_tolerance,lon_tolerance]) + z.lon.values
    lat_lim = np.array([-lat_tolerance,lat_tolerance]) + z.lat.values

    c = da.where(lon_lim[0]<da.lon,drop=True)
    c = c.where(c.lon<lon_lim[1],drop=True)

    c = c.where(lat_lim[0]<c.lat,drop=True)
    c = c.where(c.lat<lat_lim[1],drop=True)

    return c
