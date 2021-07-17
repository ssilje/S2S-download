import pandas as pd
import xarray as xr
import xskillscore as xs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from S2S.data_handler import ERA5
from S2S.process import Hindcast, Observations, Forecast, Observations_hcfc

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models
import S2S.xarray_helpers    as xh
from S2S.scoring import uncentered_acc, centered_acc
from S2S.local_configuration import config




bounds          = (0,28,55,75)

var             = 't2m'
clabel          = 'K'
t_start         = (2018,3,1)
t_end           = (2018,9,24)


clim_t_start    = (1999,1,1)
clim_t_end      = (2019,12,31)


high_res        = False

steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')
#steps           = pd.to_timedelta(np.linspace(1,46,46),'D')

print('\tProcess Forecast')
#grid_forecast = Forecast(
#                        var,
#                        t_start,
#                        t_end,
#                        bounds,
#                        high_res=high_res,
#                        steps=steps,
#                        process=False,
#                        download=False,
#                        split_work=True
#                    )

grid_forecast = Forecast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        steps=steps,
                        process=False,
                        download=False,
                        split_work=False
                    )



print('\tProcess Hindcast')
#grid_hindcast = Hindcast(
#                        var,
#                        t_start,
#                        t_end,
#                        bounds,
#                        high_res=high_res,
#                        steps=steps,
#                        process=False,
#                        download=False,
#                        split_work=True
#                    )

grid_hindcast = Hindcast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        steps=steps,  
                        process=False,
                        download=False,
                        split_work=False
                    )

## Sjekk - dataene blir berre nan. loopen under ser ok ut. 
hc_fc = []
hc_fc.append(grid_hindcast.data)
hc_fc.append(grid_forecast.data)

hindcast_full = xr.concat(hc_fc,dim='time') ## stacking the data along month dimension
hindcast_full = hindcast_full.rename(var)

era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]
grid_observations = Observations_hcfc(
                            name='Era',
                            observations=era,
                            forecast=hindcast_full,
                            process=False
                            )



reanalysis         = xh.assign_validation_time(grid_observations.data)


hindcast           = xh.assign_validation_time(grid_hindcast.data)


forecast           = xh.assign_validation_time(grid_forecast.data)

#fc_member          = list(forecast.groupby('member'))


#for n,(nn,fcdata_m) in enumerate(fc_member): # loop through each member
fcc_step = []  
hcc_step = []
re_step  = []
for lt in steps:
    fc_steps          = forecast.sel(step=pd.Timedelta(lt,'D')) #loop through each month
    hc_steps          = hindcast.sel(step=pd.Timedelta(lt,'D'))
    re_steps          = reanalysis.sel(step=pd.Timedelta(lt,'D'))   
        
    dim               = 'validation_time.month'
    
    fc_group          = list(fc_steps.groupby(dim)) 
    hc_group          = list(hc_steps.groupby(dim))
    re_group          = list(re_steps.groupby(dim))
    
    fcc_month = []
    hcc_month = []
    re_month  = [] 
    
    for m,(mf,fcdata) in enumerate(fc_group): #loop through each month
        mh,hcdata          = hc_group[m]
        mr,redata          = re_group[m]
            
        dim                 = 'validation_time.day'
        
        fc_group_day        = list(fcdata.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
        hc_group_day        = list(hcdata.groupby(dim))
        re_group_day        = list(redata.groupby(dim))
          
        fcc_day = []
        hcc_day = []
        rcc_day = []
        
        for mm,(mmf,fcdata_m) in enumerate(fc_group_day): #loop through each month
            mmh,hcdata_m          = hc_group_day[mm]
            mmr,redata_m          = re_group_day[mm]
             
            fcc = []
            for ensm in range(0,51,1):
                fc_anom = fcdata_m.sel(member=ensm) - hcdata_m.mean('time').mean('member')  
                  
                fcc.append(fc_anom)
            fcc_member = xr.concat(fcc,dim='member')    
            fcc_day.append(fcc_member)
            
            hcc = []
            for ensm in range(0,11,1):
                hc_anom = hcdata_m.mean('time').sel(member=ensm) - hcdata_m.mean('time').mean('member')  # sjekk om man bør ta mean over ens istaden og beholde kvart år
                hc_anom = hc_anom.assign_coords(time=fcc_member.time[0].values)  
                hc_anom = hcc.append(hc_anom)
            hcc_member = xr.concat(hcc,dim='member')    
            hcc_day.append(hcc_member)
            
            re_anom = redata_m - redata_m.mean('time')  
            rcc_day.append(re_anom)
            
        fcc_member_day = xr.concat(fcc_day,dim='time')
        fcc_month.append(fcc_member_day)
        
        hcc_member_day = xr.concat(hcc_day,dim='time')
        hcc_month.append(hcc_member_day)
        
        rcc_member_day = xr.concat(rcc_day,dim='time')
        re_month.append(rcc_member_day)
      
    fcc_member_day_month = xr.concat(fcc_month,dim='time')
    fcc_step.append(fcc_member_day_month)
    
    hcc_member_day_month = xr.concat(hcc_month,dim='time')
    hcc_step.append(hcc_member_day_month)
    
    re_member_day_month = xr.concat(re_month,dim='time')
    re_step.append(re_member_day_month)
    
    
    
era_anom = xr.concat(re_step,dim='step')
era_anom = era_anom.sel(time='2018')
forecast_anom = xr.concat(fcc_step,dim='step')
hindcast_anom = xr.concat(hcc_step,dim='step') 
hindcast_anom = xh.assign_validation_time(hindcast_anom)




forecast_anom_sel = forecast_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='14 days')
forecast_anom_sel = forecast_anom_sel.mean('lon').mean('lat').drop('step')
forecast_anom_df = forecast_anom_sel.to_dataframe().reset_index(level =1,drop=True).reset_index(level=0)

hindcast_anom_sel = hindcast_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='14 days')
hindcast_anom_sel = hindcast_anom_sel.mean('lon').mean('lat').drop('step')
hindcast_anom_df = hindcast_anom_sel.to_dataframe().reset_index(level =0,drop=True).reset_index(level=0)

pieces = {"fc": forecast_anom_df, "hc": hindcast_anom_df}
result = pd.concat(pieces)
test = result.reset_index(level=0)

plt.close()
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
ax = sns.boxplot(x="validation_time", y="t2m", data=test,hue='level_0', ax=axes)
x_dates = test['validation_time'].dt.strftime('%m-%d').sort_values().unique()
ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('fc_hc_test_step14.png',dpi='figure',bbox_inches='tight')
####

forecast_anom_sel = forecast_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='7 days')
forecast_anom_sel = forecast_anom_sel.mean('lon').mean('lat').drop('step')
forecast_anom_df = forecast_anom_sel.to_dataframe().reset_index(level =1,drop=True).reset_index(level=0)

hindcast_anom_sel = hindcast_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='7 days')
hindcast_anom_sel = hindcast_anom_sel.mean('lon').mean('lat').drop('step')
hindcast_anom_df = hindcast_anom_sel.to_dataframe().reset_index(level =0,drop=True).reset_index(level=0)

pieces = {"fc": forecast_anom_df, "hc": hindcast_anom_df}
result = pd.concat(pieces)
test = result.reset_index(level=0)

plt.close()
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
ax = sns.boxplot(x="validation_time", y="t2m", data=test,hue='level_0', ax=axes)
x_dates = test['validation_time'].dt.strftime('%m-%d').sort_values().unique()
ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('fc_hc_test_step7.png',dpi='figure',bbox_inches='tight')


####

forecast_anom_sel = forecast_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='21 days')
forecast_anom_sel = forecast_anom_sel.mean('lon').mean('lat').drop('step')
forecast_anom_df = forecast_anom_sel.to_dataframe().reset_index(level =1,drop=True).reset_index(level=0)

hindcast_anom_sel = hindcast_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='21 days')
hindcast_anom_sel = hindcast_anom_sel.mean('lon').mean('lat').drop('step')
hindcast_anom_df = hindcast_anom_sel.to_dataframe().reset_index(level =0,drop=True).reset_index(level=0)

pieces = {"fc": forecast_anom_df, "hc": hindcast_anom_df}
result = pd.concat(pieces)
test = result.reset_index(level=0)

plt.close()
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
ax = sns.boxplot(x="validation_time", y="t2m", data=test,hue='level_0', ax=axes)
x_dates = test['validation_time'].dt.strftime('%m-%d').sort_values().unique()
ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('fc_hc_test_step21.png',dpi='figure',bbox_inches='tight')






































anom_df = forecast_anom_df.rename(columns={'t2m':'fc-t2m'})
anom_df = hindcast_anom_df.rename(columns={'t2m':'hc-t2m'})

#pd.concat([anom_df,hindcast_anom_df.rename(columns={'t2m':'hc_t2m'})])
##anom_df.append(hc_anom,sort=False)

#pieces = {"fc": anom_df, "hc": hc_anom}
#result = pd.concat(pieces)


ax = sns.boxplot(x="validation_time", y="t2m", data=forecast_anom_df, linewidth=2.5)
x_dates = forecast_anom_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('fc_step14.png',dpi='figure',bbox_inches='tight')



ax = sns.boxplot(x="validation_time", y="t2m", data=hindcast_anom_df, linewidth=2.5)
x_dates = hindcast_anom_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('hc_step14.png',dpi='figure',bbox_inches='tight')


era_anom_sel = era_anom.sel(lat=slice(55,65), lon=slice(5,10)).sel(step='14 days')
era_anom_sel = era_anom_sel.mean('lon').mean('lat').drop('step')
era_anom_df = era_anom_sel.to_dataframe().reset_index(level =1,drop=True)

ax = sns.boxplot(x="validation_time", y="t2m", data=era_anom_df, linewidth=2.5)
x_dates = era_anom_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('era_step14.png',dpi='figure',bbox_inches='tight')


#ax = sns.boxplot(x="validation_time", y="t2m", hue="time", data=tips, linewidth=2.5) --> lag en df der hue er: hc, fc, era
.drop(
            columns=['number']
    ).reset_index(
            level=1, drop=True
    ) 

era_sel = era_anom.sel(lat="60", lon="5", method='nearest').sel(time='2018-04').sel(step='14 days').drop('time').to_dataframe()
fig, ax = plt.subplots(figsize = (12,6))

fig = sns.scatterplot(data=era_sel,
                      x="validation_time",
                      y="t2m", 
                      ax = ax)

x_dates = era_sel['validation_time'].dt.strftime('%m-%d').sort_values().unique()

ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
plt.savefig('test.png',dpi='figure',bbox_inches='tight')






## 
era_sel = era_anom.sel(lat="60", lon="5", method='nearest').sel(step='14 days')
era_sel.plot(
x=validation_time)
era_sel.plot(
x='validation_time')
plt.savefig('test.png',dpi='figure',bbox_inches='tight')
era_sel = era_anom.sel(lat="60", lon="5", method='nearest').sel(step='14 days').sel(validation_time = 2018)
era_sel = era_anom.sel(lat="60", lon="5", method='nearest').sel(step='14 days').sel(time = 2018)
era_sel = era_anom.sel(lat="60", lon="5", method='nearest').sel(step='14 days').sel(time = '2018')
era_sel
era_sel.plot(
x='validation_time')
plt.savefig('test.png',dpi='figure',bbox_inches='tight')
 plt.close()
era_sel.plot(
x='validation_time')
plt.savefig('test.png',dpi='figure',bbox_inches='tight')



      
                 # fc_anom=fc_anom.assign_coords(time_month=xlabel)
          

    
    
    fc_group_m_t          = list(fcdata_m.groupby(dim)) #loop through each month
    for nm,(nnm,fcdata_m_t) in enumerate(fc_group_m_t):
        print(nm)
        print(fcdata_m_t.dims)
  
  forecast_a        = fcdata - hindcast_mean.sel(time='2017')


## loop through each member to calculate the anomaly
  
#forecast_a        = xh.assign_validation_time(grid_forecast.data_a)
#forecast_mean     = xh.assign_validation_time(grid_forecast.mean)
#forecast_std      = xh.assign_validation_time(grid_forecast.std)



print('\tLoop through all steps')

for lt in steps:

    fc_a          = forecast_a.sel(step=pd.Timedelta(lt,'D'))
    era_a         = stacked_era_a.sel(step=pd.Timedelta(lt,'D'))
    hc_c          = hindcast_a.sel(step=pd.Timedelta(lt,'D'))
    
    
   
  
    fc_group      = list(fc_a.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    era_group     = list(era_a.groupby(dim))
    hc_group      = list(hc_c.groupby(dim))
    
