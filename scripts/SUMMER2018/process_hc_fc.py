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


for lt in steps:
      
   # fc_group_month        = list(forecast_anom_sel.groupby(dim))
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
        
        # mean over an area of norway
        fcdata_sel = fcdata.sel(lat=slice(58,65), lon=slice(5,10))
        hcdata_sel = hcdata.sel(lat=slice(58,65), lon=slice(5,10))
        eradata_sel = eradata.sel(lat=slice(58,65), lon=slice(5,10))
        
        fcdata_sel = fcdata_sel.mean('lon').mean('lat')
        hcdata_sel = hcdata_sel.mean('lon').mean('lat')
        eradata_sel = eradata_sel.mean('lon').mean('lat')
       
        fcdata_sel_df = fcdata_sel.drop('step').to_dataframe().reset_index(level = 1,drop=True).reset_index(level=0)
        hcdata_sel_df = hcdata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True).reset_index(level=0).reset_index(level=0)
        eradata_sel_df = eradata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True)
        
        
        pieces = {"fc": fcdata_sel_df, "hc": hcdata_sel_df, "era": eradata_sel_df}
        result = pd.concat(pieces)
        plot_df = result.reset_index(level=0).rename(columns={'level_0':'product'})
        
        
        plt.close()
        
        my_pal = {"fc": "g", "hc": "b", "era":"k"}
        
        
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
        ax = sns.boxplot(x="validation_time", y="t2m", data=plot_df,hue='product', ax=axes, palette=my_pal)
        x_dates = plot_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
        ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
        ax.set_ylim([-8, 8]) 
        figname = 'HC_FC_step_' + str(lt.days) + '_month_' + str(mf) + '.png'
        plt.savefig(figname,dpi='figure',bbox_inches='tight')

        
# Different plotting style
#x="validation_time", y="t2m", data=plot_df,hue='product', ax=axes, palette=my_pal)

#mm, hcdata        = hc_group[6]
#mm, eradata       = era_group[6]
#mm, fcdata        = fc_group[6]

#fcdata_sel_df = fcdata_sel.drop('step').to_dataframe().reset_index(level = 1,drop=True).reset_index(level=0)
#hcdata_sel_df = hcdata_sel.drop('step').drop('year').to_dataframe().reset_index(level=0,drop=True).reset_index(level=1,drop=True).reset_index(level=0)
#eradata_sel_df = eradata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True)

plt.close()
fig,ax=plt.subplots()
sns.boxplot(x="validation_time", y="t2m",data=hcdata_sel_df, color='b', boxprops=dict(alpha=.1),ax=ax)
sns.boxplot(x="validation_time", y="t2m",data=fcdata_sel_df, color='b',ax=ax)
ax.plot(eradata_sel_df.t2m,'r-o',linewidth=4)
plt.savefig('test_lineplot_hc_fc.png',dpi='figure',bbox_inches='tight')




#plt.close()
#fig,ax2=plt.subplots()
##sns.lineplot(x="validation_time", y="t2m",data=hcdata_sel_df,ax=ax2,color='b', boxprops=dict(alpha=.1),ax=ax2)
#sns.lineplot(x="validation_time", y="t2m",data=hcdata_sel_df,ax=ax2,color='b', alpha=.1,err_style="band",ci=100)
##sns.lineplot(x="validation_time", y="t2m",data=fcdata_sel_df,ax=ax2, color='b',err_style="bars",ci=100)
#sns.lineplot(x="validation_time", y="t2m",data=fcdata_sel_df,ax=ax2, color='b',err_style="bars")
##sns.lineplot(x="validation_time", y="t2m",data=eradata_sel_df, color='r', markers=True, dashes=True,ax=ax2)
#ax2.plot(eradata_sel_df.validation_time, eradata_sel_df.t2m,color='red', marker='o',linewidth=2,linestyle='dashed',)
#plt.savefig('test_lineplot_hc_fc_new.png',dpi='figure',bbox_inches='tight')
#bounds          = (0,28,55,75)
region = {
    'Norway': {
        'minlat': '55',  
        'maxlat': '70',
        'minlon': '5',  
        'maxlon': '20'
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
       
            fcdata_sel_df = fcdata_sel.drop('step').to_dataframe().reset_index(level = 1,drop=True).reset_index(level=0)
            hcdata_sel_df = hcdata_sel.drop('step').drop('year').to_dataframe().reset_index(level=0,drop=True).reset_index(level=1,drop=True).reset_index(level=0)
            #hcdata_sel_df = hcdata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True).reset_index(level=0).reset_index(level=0)
            eradata_sel_df = eradata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True)
        
        
            plt.close()
            fig,ax2=plt.subplots()
            sns.lineplot(x="validation_time", y="t2m",data=hcdata_sel_df,ax=ax2,color='b', alpha=.1,err_style="band",ci=100)
            sns.lineplot(x="validation_time", y="t2m",data=fcdata_sel_df,ax=ax2, color='b',err_style="bars")
            ax2.plot(eradata_sel_df.validation_time, eradata_sel_df.t2m,color='red', marker='o',linewidth=2,linestyle='dashed',)
            x_dates = eradata_sel_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
            ax2.set_xticklabels(labels=x_dates, rotation=45, ha='right')
            ax2.set_ylim([-4.5, 4.5]) 
            figname = 'HC_FC_step_' + str(lt.days) + '_month_' + str(mf) + '_' + reg + '.png'
            plt.savefig(figname,dpi='figure',bbox_inches='tight')



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
            #sns.boxplot(x="validation_time", y="t2m",data=hcdata_sel_df, color='b', boxprops=dict(alpha=.1),ax=ax)
            #sns.lineplot(x="validation_time", y="t2m",data=hcdata_sel_df,ax=ax2,color='b', alpha=.1,err_style="band",ci=100)
            ax.fill_between(eradata_sel_df.validation_time, hcmax_plot.squeeze(),hcmin_plot.squeeze(),alpha=0.1,zorder=30)
            ax.fill_between(eradata_sel_df.validation_time, fcmax_plot.squeeze(),fcmin_plot.squeeze(),alpha=0.1,zorder=30)
            sns.boxplot(x="validation_time", y="t2m",data=fcdata_sel_df, color='b',ax=ax)
            ax.plot(eradata_sel_df.t2m,'r-o',linewidth=4)
            #ax.plot(eradata_sel_df.validation_time, eradata_sel_df.t2m,color='red', marker='o',linewidth=2,linestyle='dashed',)
            x_dates = eradata_sel_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
            ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
            ax.set_ylim([-4.5, 6]) 
            figname = 'HC_FC_step_' + str(lt.days) + '_month_' + str(mf) + '_' + reg + '.png'
            plt.savefig(figname,dpi='figure',bbox_inches='tight')
            
            #plt.close()





          
            
            
            
            
            
 plt.close()
 fig,ax2=plt.subplots()
 #sns.lineplot(x="validation_time", y="t2m",data=hcdata_sel_df,ax=ax2,color='b', alpha=.1,err_style="band",ci=100)
 ax.fill_between(fcdata_sel_df.validation_time, fcdata_sel_df.t2m.max,fcdata_sel_df.t2m.min,alpha=0.1,zorder=30)
 sns.lineplot(x="validation_time", y="t2m",data=fcdata_sel_df,ax=ax2, color='b',err_style="bars")
 ax2.plot(eradata_sel_df.validation_time, eradata_sel_df.t2m,color='red', marker='o',linewidth=2,linestyle='dashed',)
 x_dates = eradata_sel_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
 ax2.set_xticklabels(labels=x_dates, rotation=45, ha='right')
 ax2.set_ylim([-4.5, 4.5]) 
 figname = 'HC_FC_step_' + str(lt.days) + '_month_' + str(mf) + '_' + reg + '.png'
 plt.savefig(figname,dpi='figure',bbox_inches='tight')       
 


        
for lt in steps:
    fc_steps          = forecast_anom.sel(step=pd.Timedelta(lt,'D')) #loop through each month  
    era_steps          = era_anom.sel(step=pd.Timedelta(lt,'D'))          
    
    dim               = 'validation_time.month'    
    
    fc_group          = list(fc_steps.groupby(dim)) 
    era_group          = list(era_steps.groupby(dim))
    
    for m,(mf,fcdata) in enumerate(fc_group): #loop through each month 
        
        me,eradata          = era_group[m]   
        
        dim                 = 'validation_time.day'    
        
        fc_group_day        = list(fcdata.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
        re_group_day        = list(eradata.groupby(dim))
          
        fcc_day = []
        rcc_day = []
        
        for mm,(mmf,fcdata_m) in enumerate(fc_group_day): #loop through each month
            mmr,redata_m          = re_group_day[mm]
            
            da = fcdata_m.mean('member')
            try:
              da = da.isel(member=0)
            except ValueError:
              pass
            
            try:
              da = da.isel(time=0)
            except ValueError:
              pass
            try:
              da = da.isel(step=0)
            except ValueError:
              pass
            
            
            p = da.transpose('lat','lon').plot(
              subplot_kws=dict(projection=ccrs.PlateCarree(),
                               facecolor="white"),
              transform=ccrs.PlateCarree(),
            )
            p.axes.coastlines()
    
            
            plt.close()
            im = fcdata_mean.drop('time').drop('step').plot( 
              x='lon',
              y='lat',
              levels=[5, 2.5, 0, 2.5, 5],
              transform=ccrs.PlateCarree(),
              cbar_kwargs={'label': 'C',
                           'ticks': [5, 2.5, 0, 2.5, 5]}
            )
            plt.savefig('test.png',dpi='figure',bbox_inches='tight')
            
            for i,ax in enumerate(im.axes.flat):
                ax.coastlines(resolution='10m', color='black',\
                      linewidth=0.2)
         
        # plotting maps
    
            #MAE_map = xs.mae(fcdata_valtime_err.mean('member',skipna=True), eradata_valtime_err, dim=[])
    
    im = varplot.plot( 
            x='lon',
            y='lat',
            col='step',
            col_wrap=3,
            levels=levels_plot,
            transform=ccrs.PlateCarree(),
            cbar_kwargs={'label': label_text,
                 'ticks': levels_cbar}
            )
  
    for i,ax in enumerate(im.axes.flat):
        ax.coastlines(resolution='10m', color='black',\
                      linewidth=0.2)
        ax.set_title(month(i))
    plt.suptitle(plot_title)
    plt.savefig(plot_save)    
            
            
            
            
            
        # mean over an area of norway
        fcdata_sel = fcdata.sel(lat=slice(58,65), lon=slice(5,10))
        hcdata_sel = hcdata.sel(lat=slice(58,65), lon=slice(5,10))
        eradata_sel = eradata.sel(lat=slice(58,65), lon=slice(5,10))
        
        fcdata_sel = fcdata_sel.mean('lon').mean('lat')
        hcdata_sel = hcdata_sel.mean('lon').mean('lat')
        eradata_sel = eradata_sel.mean('lon').mean('lat')
       
        fcdata_sel_df = fcdata_sel.drop('step').to_dataframe().reset_index(level = 1,drop=True).reset_index(level=0)
        hcdata_sel_df = hcdata_sel.drop('step').drop('year').to_dataframe().reset_index(level=0,drop=True).reset_index(level=1,drop=True).reset_index(level=0)
        #hcdata_sel_df = hcdata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True).reset_index(level=0).reset_index(level=0)
        eradata_sel_df = eradata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True)
        

        
    def plot_months(
        varplot,
        levels_plot,
        label_text,
        levels_cbar,
        plot_title,
        plot_save,
):
    # Sjekk her: plottet blir det samme om eg brukar transpose eller ikkje. 
    #im = skill_score_at_lt.skill.transpose('lon','lat','time_month').plot(
    im = varplot.plot( 
        x='lon',
        y='lat',
        col='time_month',
        col_wrap=3,
        levels=levels_plot,
        transform=ccrs.PlateCarree(),
        cbar_kwargs={'label': label_text,
                 'ticks': levels_cbar}
    )
  
    for i,ax in enumerate(im.axes.flat):
        ax.coastlines(resolution='10m', color='black',\
                      linewidth=0.2)
        ax.set_title(month(i))
    plt.suptitle(plot_title)
    plt.savefig(plot_save)    







