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
    

forecast_length = pd.date_range(start="20200505", periods=46) # forecasts start Monday

f, ax = plt.subplots(1, 1)
#for ens in range(1,11):
#    print(ens)
    #ax.plot(dataopen_pf_hc.valid_time[dataopen_pf_hc.valid_time.index.get_level_values('number') == ens], dataopen_pf_hc.sav300[dataopen_pf_hc.sav300.index.get_level_values('number') == ens],color='b',linewidth=0.2)
#    ax.plot(dataopen_pf_hc.index.get_level_values('step'), dataopen_pf_hc.sav300[dataopen_pf_hc.sav300.index.get_level_values('number') == ens],color='b',linewidth=0.2)

#ax.plot(dataopen_cf.valid_time, dataopen_cf_hc.sav300, linewidth=0.2, color='b')


for ens in range(1,51):
    print(ens)
    #ax.plot(dataopen_pf.valid_time[dataopen_pf.valid_time.index.get_level_values('number') == ens], dataopen_pf.sav300[dataopen_pf.sav300.index.get_level_values('number') == ens],color='r',linewidth=0.5)
    ax.plot(dataopen_pf.valid_time[dataopen_pf.valid_time.index.get_level_values('number') == ens], dataopen_pf.sav300[dataopen_pf.sav300.index.get_level_values('number') == ens] ,color='r',linewidth=0.5)
   # ax.plot(dataopen_pf.valid_time[dataopen_pf.valid_time.index.get_level_values('number') == ens].index.get_level_values('step'), dataopen_pf.sav300[dataopen_pf.sav300.index.get_level_values('number') == ens] ,color='r',linewidth=0.5)
#ax.plot(dataopen_cf.valid_time, dataopen_cf.sav300, linewidth=2, color='k')
ax.plot(dataopen_cf.index.get_level_values('step'), dataopen_cf.sav300, linewidth=2, color='k')

xfmt = mdates.DateFormatter('%m%d')
ax.xaxis.set_major_formatter(xfmt)

ax.set_xlabel('time')
ax.set_ylabel('SALINITY')
#ax.set_title('SST July', fontsize=16)



#ERA5_BR.SST_std.plot()-ERA5_BR.SST
#ERA5_BR.groupby(ERA5_BR_df.index.1).mean().plot()
f.savefig('SALINITY.png')


#for idata_all, rdata_all in data_all.iterrows(): 
#    for iera, rera in ERA5_BR_daily.iterrows():   
#        if rdata_all.valid_time == rera.time
#           rdata_all.sst
##    if (row.valid_time==ERA5_BR_daily.index) is True:
 #   ERA5_BR_daily.index==row.valid_time
