import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from openpyxl import Workbook

var_short = 'sal' 
var = 'sav300'
cycle = 'CY46R1'
dirbase_S2S = '/nird/projects/NS9001K/sso102/S2S/DATA/grib'
lead_time = np.arange(1,47)
fcyear = 2020
fcmonth=5
fcday=4
# Bergen
lat = 65
lon = 0
dates_fcycle = pd.date_range(start='%s-%s-%s'%(fcyear,fcmonth,fcday), periods=1, freq="7D") # forecasts start Monday

def read_grib(dirbase_S2S,product,ftype,d,lat,lon):
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')
    dS2S = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,ftype,product,'.grb')
    print('reading file:')
    print(dS2S)
    dataopen = xr.open_dataset(dS2S,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    return dataopen

def read_grib_cf_pf(dirbase_S2S,product,d,lat,lon):
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')
    
    dS2S_cf = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'cf',product,'.grb')
    dS2S_pf = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'pf',product,'.grb')
    
    print('reading file:')
    print(dS2S_pf)   
    dataopen = xr.open_dataset(dS2S_pf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point

    print('reading file:')
    print(dS2S_cf)   
    dataopen_cf = xr.open_dataset(dS2S_cf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    newdata = pd.Series(
        {"time": [dataopen_cf.time], 
        "oceanLayer": [dataopen_cf.oceanLayer], 
        "latitude": [dataopen_cf.latitude], 
        "longitude": [dataopen_cf.longitude], 
        "valid_time": [dataopen_cf.valid_time], 
        "sav300": [dataopen_cf.sav300]},
        index=pd.MultiIndex.from_product([["number", "step"], [dataopen_cf.number,dataopen_cf.index.get_level_values('step')]]),
   
    
    #test_dataopen = dataopen.join(dataopen_cf.set_index('number'), on='number')
    index = pd.MultiIndex.from_tuples(dataopen_cf['number'], names=["step","number"])
    data_stats = pd.DataFrame({"mean": [data_mean], "std": [data_std]}, index=[step])
    data_stats.index.name = 'step'    
    
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

def calc_stats_lead_time_cf_pf(dataopen_cf,dataopen_pf,step,var,ftype):
    data_pf_cf = pd.concat([dataopen_cf_pf.loc[(step,slice(None)),var], dataopen_cf_pf.loc[(slice(None),step,slice(None)),var])
    
    if ftype == 'cf':
        data_mean = dataopen.loc[(step,slice(None)),var].mean()
        data_std = dataopen.loc[(step,slice(None)),var].std()
    if ftype == 'pf':
        data_mean = dataopen.loc[(slice(None),step,slice(None)),var].mean()
        data_std = dataopen.loc[(slice(None),step,slice(None)),var].std()
        
    data_stats = pd.DataFrame({"mean": [data_mean], "std": [data_std]}, index=[step])
    data_stats.index.name = 'step'
    
    return data_stats
    

for idate in dates_fcycle: 

    d = idate.strftime('%Y-%m-%d')
    dataopen_cf = read_grib(dirbase_S2S,'forecast','cf',d,lat,lon) #product, ftype, lat, lon
    dataopen_pf = read_grib(dirbase_S2S,'forecast','pf',d,lat,lon) #product, ftype, lat, lon
    dataopen_cf_hc = read_grib(dirbase_S2S,'hindcast','cf',d,lat,lon) #product, ftype, lat, lon
    dataopen_pf_hc = read_grib(dirbase_S2S,'hindcast','pf',d,lat,lon) #product, ftype, lat, lon
   
#ns = set(dataopen_cf_hc.index.get_level_values('step').days) # set getst the unique values

# Pandas timedeltas: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.TimedeltaIndex.html

for nstep in set(dataopen_cf_hc.index.get_level_values('step').days): # set getst the unique values
    #step = '%s%s'%(nstep,' days')
    print(nstep)
    print('%s%s'%(nstep,' days'))
    if nstep == 1: 
        data_stats_cf_hc = calc_stats_lead_time(dataopen_cf_hc,'%s%s'%(nstep,' days'),var,'cf')
        data_stats_pf_hc = calc_stats_lead_time(dataopen_pf_hc,'%s%s'%(nstep,' days'),var,'pf')
 
    else:
        data_stats_cf_hc = pd.concat([data_stats_cf_hc, calc_stats_lead_time(dataopen_cf_hc,'%s%s'%(nstep,' days'),var,'cf')])
        data_stats_pf_hc = pd.concat([data_stats_pf_hc, calc_stats_lead_time(dataopen_pf_hc,'%s%s'%(nstep,' days'),var,'pf')])

        
f, axes = plt.subplots(1, 1)
ax = axes

mmm = data_stats_cf_hc.mean

# this loops through each row in the array
for y in tas_anom:
    # plot here
    pass

h = ax.plot(year, mmm, lw=2)

ax.set_ylabel('Tanom [°C]')
ax.set_xlabel('Time')

ax.set_title('Global mean temperature', fontsize=14)

ax.axhline(0, color='0.1', lw=0.5)

ax.axvspan(1971, 2000, color='0.75')

plt.show()


        
        

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

