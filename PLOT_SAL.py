import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from openpyxl import Workbook

var_short = 'sal' 

ftype = 'pf'


cycle = 'CY46R1'

# Bergen
lat = 65
lon = 0

dirbase_S2S = '/nird/projects/NS9001K/sso102/S2S/DATA/grib'

dates_monday = pd.date_range("20200504", periods=1, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20190704", periods=1, freq="7D") # forecasts start Thursday

dates_fcycle = dates_monday

for idate in dates_fcycle: 
#for idate in dates_monday: 
    d = idate.strftime('%Y-%m-%d')
    product = 'forecast' # forecast
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')
    dS2S_cf = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'cf',product,'.grb')
    dataopen_cf = xr.open_dataset(dS2S_cf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    #dataopen_cf.reset_index(inplace=True)  
    dS2S_pf = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'pf',product,'.grb')
    dataopen_pf = xr.open_dataset(dS2S_pf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    #dataopen_pf.reset_index(inplace=True)  
    
    product = 'hindcast' # forecast
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')
    dS2S_cf_hc = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'cf',product,'.grb')
    dataopen_cf_hc = xr.open_dataset(dS2S_cf_hc,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    #dataopen_cf_hc.reset_index(inplace=True)  
    dS2S_pf_hc = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'pf',product,'.grb')
    dataopen_pf_hc = xr.open_dataset(dS2S_pf_hc,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    #dataopen_pf_hc.reset_index(inplace=True) 
    


lead_time = np.arange(1,47)
f, ax = plt.subplots(1, 1)
for y in range(2000,2020):
    year = '%s'%(y)
    d = pd.date_range(start='%s-%s-%s'%(y,'05','04'), periods=46) # forecasts start Monday

    for ens in range(1,11):    
        ax.plot(lead_time, dataopen_pf_hc.sav300.loc[(ens,slice(None), year)],color='b',linewidth=0.2)
    ax.plot(lead_time, dataopen_cf_hc.sav300.loc[((slice(None),d.strftime('%Y%m%d')))],color='k',linewidth=0.2)

for ens in range(1,51):
    ax.plot(lead_time, dataopen_pf.sav300.loc[(ens)] ,color='r',linewidth=0.5)
ax.plot(lead_time, dataopen_cf.sav300 ,color='k',linewidth=1)


ax.set_xlabel('lead time')
ax.set_ylabel('SALINITY')



plt.show()

