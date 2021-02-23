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

dates_monday = pd.date_range("20190701", periods=1, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20190704", periods=1, freq="7D") # forecasts start Thursday

dates_fcycle = dates_monday.union(dates_thursday)   


for idate in dates_fcycle: 
#for idate in dates_monday: 
    d = idate.strftime('%Y-%m-%d')
    
    dS2S_cf = '%s/%s/%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'cf','.grb')
    dataopen_cf = xr.open_dataset(dS2S_cf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    dS2S_pf = '%s/%s/%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'pf','.grb')
    dataopen_pf = xr.open_dataset(dS2S_pf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    #dataopen_pf.reset_index(inplace=True)
    #dataopen_cf.reset_index(inplace=True)  
    data_append = dataopen_cf.reset_index(inplace=True).append(dataopen_pf).reset_index(inplace=True)  
    
  #  data_append.reset_index(inplace=True)  
    if d == dates_monday[0].strftime('%Y-%m-%d'):
        
        data_all = data_append
    else:                             
        data_all = data_all.append(data_append) 
    
    #data_all.reset_index(inplace=True)  
    print('data_all')
    print(data_all.head(50)) # print the 20 first lines
    
    
    
