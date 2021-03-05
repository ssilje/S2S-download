
import xarray as xr
#import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime

import globals as stosgl


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

for idate in dates_fcycle: 

    d = idate.strftime('%Y-%m-%d')
    dataopen_fc =stosgl.read_grib_cf_pf(dirbase_S2S,'forecast',d,lat,lon,var_short,cycle)
    dataopen_hc = stosgl.read_grib_cf_pf(dirbase_S2S,'hindcast',d,lat,lon,var_short,cycle)
    

for nstep in set(dataopen_hc.index.get_level_values('step').days): # set gets the unique values
    print('%s%s'%(nstep,' days'))
    if nstep == 1: 
     
        data_stats_hc = stosgl.calc_stats_lead_time_cf_pf(dataopen_hc,'%s%s'%(nstep,' days'),var)
        data_stats_fc = stosgl.calc_stats_lead_time_cf_pf(dataopen_fc,'%s%s'%(nstep,' days'),var)
 
    else:
        data_stats_hc = pd.concat([data_stats_hc, stosgl.calc_stats_lead_time_cf_pf(dataopen_hc,'%s%s'%(nstep,' days'),var)])
        data_stats_fc = pd.concat([data_stats_fc, stosgl.calc_stats_lead_time_cf_pf(dataopen_fc,'%s%s'%(nstep,' days'),var)])
        

f, ax = plt.subplots(1, 1)

ax.plot(lead_time, data_stats_hc.data_mean.loc[(ens,slice(None), year)],color='gray',linewidth=0.2)
ax.plot(lead_time, dataopen_cf_hc.sav300.loc[((slice(None),d.strftime('%Y%m%d')))],color='gray',linewidth=0.2)

for ens in range(1,51):
    ax.plot(lead_time, dataopen_pf.sav300.loc[(ens)] ,color='blue',linewidth=0.5)
ax.plot(lead_time, dataopen_cf.sav300 ,color='k',linewidth=1)
ax.plot(lead_time, dataopen_pf.sav300.loc[(ens)] ,color='blue',linewidth=0.5)

ax.set_xlabel('lead time')
ax.set_ylabel('SALINITY')

plt.show()
