import xarray as xr
import pandas as pd

# # get observations
# path     = '/nird/projects/NS9853K/DATA/norkyst800/'
# filename = 'norkyst800_sst_*_daily_mean_at-BW-loc.nc'
#
# ds1 = xr.open_dataset(
#                 path + 'norkyst800_sst_2016-11-12_daily_mean_at-BW-loc.nc'
#             )
# print(ds1)
# exit()
# ds = xr.open_mfdataset( path + filename )
#
# print(ds)

##################################################################
start_date = '2012-06-27'
end_date   = '2019-02-26'

date_range = pd.date_range(start=start_date,end=end_date,freq='D')

for n,date in enumerate(date_range):

    filename_nird = 'norkyst800_sst_'+date.strftime('%Y-%m-%d')+\
                        '_daily_mean_at-BW-loc.nc'

    ds = xr.open_dataset( path + filename_nird )

    print(ds)
    exit()
