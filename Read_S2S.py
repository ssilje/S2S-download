import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime

var_short = 't2m' 
ftype = 'pf'
product = 'hindcast' # forecast
cycle = 'CY43R3_CY45R1'
# Bergen
lat = 60.23
lon = 5.19

#ds_grib = xr.open_dataset('data/weather/min_max_t_and_wind_gust.grib',engine='cfgrib')

dirbase_S2S = '/nird/projects/NS9853K/DATA/S2S/SUMMER2018/'
dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')

dates_monday = pd.date_range("20180426", periods=1, freq="7D") # forecasts start Monday
#dates_monday = pd.date_range("20180426", periods=20, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20180430", periods=1, freq="7D") # forecasts start Thursday
dates_fcycle = dates_monday.union(dates_thursday)  

for idate in dates_thursday: 
    d = idate.strftime('%Y-%m-%d')
    dS2S = '%s/%s/%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,ftype,'.grb')
    dataopen = xr.open_dataset(dS2S,engine='cfgrib')
    S2S_BR_daily = dataopen.sel(latitude=lat, longitude=lon, method='nearest').to_dataframe()
    print(S2S_BR_daily.head(20))
    
    dates_hc = pd.date_range((idate-pd.DateOffset(years=20)), periods=20, freq="AS-JUL") #20 years hindcast
    print(dates_hc)
    
    
