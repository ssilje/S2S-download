import netCDF4
import xarray as xr
import pandas as pd
import os
import numpy as np

from S2S.data_handler import BarentsWatch

var = 'sst'
# bw_obs = BarentsWatch().load('all',no=420).sortby('time')[var].isel(time=0)

download_full = False

url  = 'https://thredds.met.no/thredds/dodsC/sea/norkyst800mv0_24h/'

path = '/nird/projects/NS9853K/DATA/norkyst800/'

if not os.path.exists(path):
    os.mkdir(path)

start_date = '2012-06-27'
end_date   = '2019-02-26'

date_range = pd.date_range(start=start_date,end=end_date,freq='D')

fname  = 'NorKyst-800m_ZDEPTHS_avg.an.'
add_on = '12.nc'

for n,date in enumerate(date_range):

    filename = fname + date.strftime('%Y%m%d') + add_on
    filename_nird = 'norkyst800_sst_'+date.strftime('%Y-%m-%d')+'_daily_mean.nc'

    try:

        with netCDF4.Dataset(url+filename) as file:

            print(filename)

            ds = xr.open_dataset(
                                    xr.backends.NetCDF4DataStore(file),
                                    decode_times=False
                                )

            data = ds.temperature.sel(depth=0.)

            if download_full:
                data.to_netcdf(path+filename_nird)

            if n==0:
                ds.sel(depth=0.).mask.to_netcdf(path+'norkyst800_landmask.nc')

    except OSError:

        pass
