import xarray as xr

# get observations
path     = '/nird/projects/NS9853K/DATA/norkyst800/'
filename = 'norkyst800_sst_*_daily_mean_at-BW-loc.nc'

ds1 = xr.open_dataset(
                path + 'norkyst800_sst_2016-11-12_daily_mean_at-BW-loc.nc'
            )
print(ds1)
exit()
ds = xr.open_mfdataset( path + filename )

print(ds)
