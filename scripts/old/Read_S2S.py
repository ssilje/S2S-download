import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime

var_short = 'sst' 
ftype = 'pf'
product = 'hindcast' # forecast

cycle = 'CY46R1'

# Bergen
lat = 60.23
lon = 5.19



dirbase_S2S = '/nird/projects/NS9001K/sso102/S2S/DATA/grib'
dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')

dates_monday = pd.date_range("20190701", periods=52, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20190704", periods=52, freq="7D") # forecasts start Thursday
dates_fcycle = dates_monday.union(dates_thursday)   


for idate in dates_fcycle: 
    d = idate.strftime('%Y-%m-%d')
    dS2S = '%s/%s/%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,ftype,'.grb')
    dataopen = xr.open_dataset(dS2S,engine='cfgrib')
    S2S_BR_daily = dataopen.sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    
    print(S2S_BR_daily.head(20)) # print the 20 first lines
