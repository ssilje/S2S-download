import xarray as xr

# get observations
path     = '/nird/projects/NS9853K/DATA/norkyst800/'
filename = 'norkyst800_sst_*_daily_mean_at-BW-loc.nc'

ds = xr.open_mfdataset( path + filename )

print(ds)
