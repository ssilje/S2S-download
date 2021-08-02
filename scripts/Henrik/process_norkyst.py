import netCDF4
import xarray as xr
import pandas as pd
import os
import numpy as np

from S2S.data_handler import BarentsWatch

def lin_weight(data,dist,r):
    """
    Linearly weighted mean where the closest points are given more weight.
    """
    weights = 1 - dist/r
    return ( data * weights ).sum() / weights.sum()

def distance(lat1, lon1, lat2, lon2):
    """
    Returns the great circle distance between lat1,lon1 and lat2,lon2.
    Function is positive definite.
    """
    p = np.pi/180
    a = 0.5 - np.cos((lat2-lat1)*p)/2 +\
        np.cos(lat1*p) * np.cos(lat2*p) * (1-np.cos((lon2-lon1)*p))/2
    return 12742 * np.arcsin(np.sqrt(a))

def norkyst_to_location(data,lat,lon,landmask,reflat,reflon):

    appr_data = []
    rs        = []

    landmask  = landmask.astype('int')

    for rlon,rlat in zip(reflon,reflat):

        dist = distance(rlat,rlon,lat,lon)

        over_land     = 0
        target_radius = 0.
        while not over_land:

            target_radius += 1.

            idx  = dist < target_radius

            lm   = landmask[idx]

            over_land = lm.sum()

        lm = lm.astype('bool')

        # wr: within radius
        data_wr = data[idx][lm]
        dist_wr = dist[idx][lm]

        appr_data.append( lin_weight( data_wr, dist_wr, target_radius ) )
        rs.append( target_radius )

    return np.array(appr_data),np.array(rs)

var = 'sst'
bw_obs = BarentsWatch().load('all',no=0).sortby('time')[var].isel(time=0)

path = '/nird/projects/NS9853K/DATA/norkyst800/'

if not os.path.exists(path):
    os.mkdir(path)

url    = 'https://thredds.met.no/thredds/dodsC/sea/norkyst800mv0_24h/'
fname  = 'NorKyst-800m_ZDEPTHS_avg.an.'
add_on = '12.nc'

start_date = '2012-06-27'
end_date   = '2019-02-26'

download = True

date_range = pd.date_range(start=start_date,end=end_date,freq='D')

landmask = xr.open_dataset(path+'norkyst800_landmask.nc').mask

big_chunk = []

for n,date in enumerate(date_range):

    filename = fname + date.strftime('%Y%m%d') + add_on

    filename_nird = 'norkyst800_sst_'+date.strftime('%Y-%m-%d')+\
                        '_daily_mean_at-BW-loc.nc'

    try:

        with netCDF4.Dataset(url+filename) as file:

            print(filename)

            ds = xr.open_dataset(
                                    xr.backends.NetCDF4DataStore(file),
                                    decode_times=False
                                )

            data     = ds.temperature.sel(depth=0.).rename(var)
            landmask = ds.sel(depth=0.).mask

            out,r = xr.apply_ufunc(
                    norkyst_to_location,
                    data,
                    data.latitude,
                    data.longitude,
                    landmask,
                    bw_obs.lat,
                    bw_obs.lon,
                    input_core_dims=[
                                        ['X','Y'],['X','Y'],['X','Y'],['X','Y'],
                                        ['location'],['location'],
                                    ],
                    output_core_dims=[['location'],['location']],
                    vectorize=True
                    # dask='parallelized'
                )

            out = out.assign_coords(radius=r)
            out = out.assign_coords(time=pd.Timestamp(date))
            print(out)
            exit()

            if download:
                out.to_netcdf(filename_nird)

    except OSError:

        pass
