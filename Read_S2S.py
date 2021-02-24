import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from openpyxl import Workbook

var_long='sea_surface_temperature' 
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
    dataopen_pf.reset_index(inplace=True)
    dataopen_cf.reset_index(inplace=True)  
    data_append = dataopen_cf.append(dataopen_pf)
    
  #  data_append.reset_index(inplace=True)  
    if d == dates_monday[0].strftime('%Y-%m-%d'):
        
        data_all = data_append
    else:                             
        data_all = data_all.append(data_append) 
    
    #data_all.reset_index(inplace=True)  
   # print('data_all')
   # print(data_all.head(50)) # print the 20 first lines
print('data_all')
print(data_all.head(50)) # print the 20 first lines

dirbase = '/nird/projects/NS9853K/DATA/SFE/ERA_daily_nc'
dates_era = pd.date_range(start='19990701', end='19991001', freq="D")
for i,d in enumerate(dates_era):
    dERA5 = '%s/%s_%s%s'%(dirbase,var_long,d.strftime('%Y%m%d'),'.nc')
    dataopen = xr.open_dataset(dERA5)
    if i == 0:
        ERA5_BR_daily = dataopen.sst.sel(lat=lat, lon=lon, method='nearest').resample(time='D').mean().to_dataframe() # Maa interpolere 
    else:
        ERA5_BR_daily = pd.concat([ERA5_BR_daily, dataopen.sst.sel(lat=lat, lon=lon, method='nearest').resample(time='D').mean().to_dataframe()])
                
                
#data_all.to_excel("output.xlsx")  
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_excel.html
    
