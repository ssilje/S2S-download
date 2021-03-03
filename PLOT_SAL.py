import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from openpyxl import Workbook

var_short = 'sal' 
cycle = 'CY46R1'
dirbase_S2S = '/nird/projects/NS9001K/sso102/S2S/DATA/grib'
lead_time = np.arange(1,47)
fcyear = 2020
fcmonth=5
fcday=4
# Bergen
lat = 65
lon = 0

def read_grib(dirbase_S2S,product,ftype,d,lat,lon):
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')
    dS2S = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,ftype,product,'.grb')
    print('reading file:')
    print(dS2S)
    dataopen = xr.open_dataset(dS2S,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    return dataopen

def calc_stats_lead_time(dataopen,step,var,ftype):
    if ftype == 'cf':
        data_mean = dataopen.loc[(step,slice(None)),var].mean()
        data_std = dataopen.loc[(step,slice(None)),var].std()
    if ftype == 'pf':
        data_mean = dataopen.loc[(slice(None),step,slice(None)),var].mean()
        data_std = dataopen.loc[(slice(None),step,slice(None)),var].std()
        
    data_stats = pd.DataFrame({"mean": [data_mean], "std": [data_std]}, index=[step])
    data_stats.index.name = 'step'
    
    return data_stats
    

dates_fcycle = pd.date_range(start='%s-%s-%s'%(fcyear,fcmonth,fcday), periods=1, freq="7D") # forecasts start Monday

for idate in dates_fcycle: 

    d = idate.strftime('%Y-%m-%d')
    dataopen_cf = read_grib(dirbase_S2S,'forecast','cf',d,lat,lon) #product, ftype, lat, lon
    dataopen_pf = read_grib(dirbase_S2S,'forecast','pf',d,lat,lon) #product, ftype, lat, lon
    dataopen_cf_hc = read_grib(dirbase_S2S,'hindcast','cf',d,lat,lon) #product, ftype, lat, lon
    dataopen_pf_hc = read_grib(dirbase_S2S,'hindcast','pf',d,lat,lon) #product, ftype, lat, lon
   
nstep = set(dataopen_cf_hc.index.get_level_values('step')) # set getst the unique values

# Pandas timedeltas: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.TimedeltaIndex.html

calc_stats_lead_time(dataopen_cf_hc,step,var,'cf')


f, ax = plt.subplots(1, 1)
syr=dataopen_cf_hc.valid_time.min().year
for y in range(syr,fcyear): # need to check for the hcast years
    year = '%s'%(y)
    d = pd.date_range(start='%s-%s-%s'%(y,fcmonth,fcday), periods=46) # forecasts start Monday

    for ens in range(1,11):    
        ax.plot(lead_time, dataopen_pf_hc.sav300.loc[(ens,slice(None), year)],color='gray',linewidth=0.2)
    ax.plot(lead_time, dataopen_cf_hc.sav300.loc[((slice(None),d.strftime('%Y%m%d')))],color='gray',linewidth=0.2)

for ens in range(1,51):
    ax.plot(lead_time, dataopen_pf.sav300.loc[(ens)] ,color='blue',linewidth=0.5)
ax.plot(lead_time, dataopen_cf.sav300 ,color='k',linewidth=1)
ax.plot(lead_time, dataopen_pf.sav300.loc[(ens)] ,color='blue',linewidth=0.5)

ax.set_xlabel('lead time')
ax.set_ylabel('SALINITY')

plt.show()

