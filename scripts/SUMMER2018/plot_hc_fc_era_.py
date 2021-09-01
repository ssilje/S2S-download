import pandas as pd
import xarray as xr
import xskillscore as xs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns

from S2S.data_handler import ERA5
from S2S.process import Hindcast, Observations, Forecast, Observations_hcfc

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models
import S2S.xarray_helpers    as xh
from S2S.scoring import uncentered_acc, centered_acc
from S2S.local_configuration import config




#bounds          = (0,28,40,75)
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
                        process=True,
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
                hcdata_m_years = list(hcdata_m.groupby('validation_time.year'))
                hcc_year = []
                for yy,(yyh,hcdata_year) in enumerate(hcdata_m_years): #loop through each year
                    
                   # hc_anom = hcdata_m.mean('time').sel(member=ensm) - hcdata_m.mean('time').mean('member')  # sjekk om man bør ta mean over ens istaden og beholde kvart år
                   # hc_anom = hc_anom.assign_coords(time=fcc_member.time[0].values)  
                   # hc_anom = hcc.append(hc_anom)
  
                    hc_anom = hcdata_year.sel(member=ensm).drop('time').drop('validation_time').squeeze() - hcdata_m.mean('time').mean('member')
                    hc_anom = hc_anom.assign_coords(time=fcc_member.time[0].values)
                    hc_anom = hc_anom.assign_coords(year=yyh)
                    hc_anom = hcc_year.append(hc_anom)

                hcc_m_year = xr.concat(hcc_year,dim='year')
                
                hc_member_anom = hcc.append(hcc_m_year)
            
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


region = {
    'Norway': {
        'minlat': '58',  
        'maxlat': '64',
        'minlon': '5',  
        'maxlon': '12'
            },
    'MEU': {
        'minlat': '45',  
        'maxlat': '55',
        'minlon': '0',  
        'maxlon': '16'
            },
}

for lt in steps:
    fc_steps          = forecast_anom.sel(step=pd.Timedelta(lt,'D')) #loop through each month
    hc_steps          = hindcast_anom.sel(step=pd.Timedelta(lt,'D'))
    era_steps          = era_anom.sel(step=pd.Timedelta(lt,'D'))         
    
    dim               = 'validation_time.month'    
    
    fc_group          = list(fc_steps.groupby(dim)) 
    hc_group          = list(hc_steps.groupby(dim))
    era_group          = list(era_steps.groupby(dim))
  
    for m,(mf,fcdata) in enumerate(fc_group): #loop through each month
        mh,hcdata          = hc_group[m]
        me,eradata          = era_group[m]
        
        for reg in (
          'Norway',
        #  'MEU'
        ):
            fcdata_sel = fcdata.sel(lat=slice(int(region[reg]['minlat']),int(region[reg]['maxlat'])), lon=slice(int(region[reg]['minlon']),int(region[reg]['maxlon'])))
            hcdata_sel = hcdata.sel(lat=slice(int(region[reg]['minlat']),int(region[reg]['maxlat'])), lon=slice(int(region[reg]['minlon']),int(region[reg]['maxlon'])))
            eradata_sel = eradata.sel(lat=slice(int(region[reg]['minlat']),int(region[reg]['maxlat'])), lon=slice(int(region[reg]['minlon']),int(region[reg]['maxlon'])))
        
            fcdata_sel = fcdata_sel.mean('lon').mean('lat')
            hcdata_sel = hcdata_sel.mean('lon').mean('lat')
            eradata_sel = eradata_sel.mean('lon').mean('lat')
            
            dim                 = 'validation_time.day'
            hc_group_day        = list(hcdata_sel.groupby(dim))
            fc_group_day        = list(fcdata_sel.groupby(dim))
            
            hcmax = []
            hcmin = []
            
            fcmax = []
            fcmin = []
            fc50 = []
            fc75 = []
            fc25 = []
            
            for hcm,(hcmm,hcdata_sel_day) in enumerate(hc_group_day): #loop through each month
                fcm, fcdata_sel_day = fc_group_day[hcm]
                
                max_hc = hcdata_sel_day.max('year').max('member').drop('time').drop('step').drop('validation_time')
                min_hc = hcdata_sel_day.min('year').min('member').drop('time').drop('step').drop('validation_time')
                hcmax.append(max_hc)
                hcmin.append(min_hc)
                
                max_fc = fcdata_sel_day.quantile(1)
                min_fc = fcdata_sel_day.quantile(0)
                p50_fc = fcdata_sel_day.quantile(0.5)
                p25_fc = fcdata_sel_day.quantile(0.25)
                p75_fc = fcdata_sel_day.quantile(0.75)
                
                fcmax.append(max_fc)
                fcmin.append(min_fc)
                fc50.append(p50_fc)
                fc75.append(p75_fc)
                fc25.append(p25_fc)
                
                
            hcmax_plot = xr.concat(hcmax,dim='validation_time') 
            hcmin_plot = xr.concat(hcmin,dim='validation_time') 
            fcmax_plot = xr.concat(fcmax,dim='validation_time')
            fcmin_plot = xr.concat(fcmin,dim='validation_time')
            fc50_plot = xr.concat(fc50,dim='validation_time')
            fc75_plot = xr.concat(fc75,dim='validation_time')
            fc25_plot = xr.concat(fc25,dim='validation_time')
                 
            fcdata_sel_df = fcdata_sel.drop('step').to_dataframe().reset_index(level = 1,drop=True).reset_index(level=0)
            hcdata_sel_df = hcdata_sel.drop('step').drop('year').to_dataframe().reset_index(level=0,drop=True).reset_index(level=1,drop=True).reset_index(level=0)
            #hcdata_sel_df = hcdata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True).reset_index(level=0).reset_index(level=0)
            eradata_sel_df = eradata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True)
        
        
            plt.close()
            fig,ax=plt.subplots()
            ax.plot(eradata_sel_df.validation_time,np.zeros(eradata_sel_df.validation_time.shape[0]),alpha=0.1)
            hclabel = ax.fill_between(eradata_sel_df.validation_time, hcmax_plot.squeeze(),hcmin_plot.squeeze(),alpha=0.1,zorder=30, facecolor='gray', label='hindcast')
            fclabel = ax.fill_between(eradata_sel_df.validation_time, fcmax_plot.squeeze(),fcmin_plot.squeeze(),alpha=0.1,zorder=30, facecolor='blue', label='forecast')
            fclabel_25_75 = ax.fill_between(eradata_sel_df.validation_time, fc75_plot.squeeze(),fc25_plot.squeeze(),alpha=0.1,zorder=30, facecolor='purple',label='forecast Q1-Q3')
            fclabel50 = ax.plot(eradata_sel_df.validation_time,fc50_plot.squeeze(),color='purple', label='forecast median')
            era = ax.plot(eradata_sel_df.validation_time,eradata_sel_df.t2m, color='red',label='ERA5')
            ax.legend(loc='lower right')
            x_dates = eradata_sel_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
            ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
            ax.set_ylim([-8, 8]) 
            figname = 'HC_FC_step_' + str(lt.days) + '_month_' + str(mf) + '_' + reg + '_full_ens.png'
            plt.savefig(figname,dpi='figure',bbox_inches='tight')
            
            
            






          
            
            
       
            
 
