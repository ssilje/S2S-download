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
        
        dim                 = 'validation_time.day'
        hc_group_day        = list(hcdata.groupby(dim))
        fc_group_day        = list(fcdata.groupby(dim))
        era_group_day       = list(eradata.groupby(dim))
            
        hcmax = []
        hcmin = []
            
        fcmax = []
        fcmin = []
        fc50 = []
        fc75 = []
        fc25 = []
        plotdata_test = []
        plotdata = []
            
        for hcm,(hcmm,hcdata_sel_day) in enumerate(hc_group_day): #loop through each month
            fcm, fcdata_sel_day = fc_group_day[hcm]
            ecm, era_sel_day = era_group_day[hcm]
            
            plotdata_test = []
            plotdata = []    
            max_hc = hcdata_sel_day.max('year').max('member').drop('time').drop('step').drop('validation_time')
            min_hc = hcdata_sel_day.min('year').min('member').drop('time').drop('step').drop('validation_time')
            hcmax.append(max_hc)
            hcmin.append(min_hc)
            
            
                
            max_fc = fcdata_sel_day.quantile(1,dim='member')
            plotdata_test.append(max_fc)
            min_fc = fcdata_sel_day.quantile(0,dim='member')
            plotdata_test.append(min_fc)
            p50_fc = fcdata_sel_day.quantile(0.5,dim='member')
            plotdata_test.append(p50_fc)
            p25_fc = fcdata_sel_day.quantile(0.25,dim='member')
            plotdata_test.append(p25_fc)
            p75_fc = fcdata_sel_day.quantile(0.75,dim='member')
            plotdata_test.append(p75_fc)
            
            era_sel_day = era_sel_day.assign_coords(quantile='era')
            plotdata_test.append(era_sel_day)
            
            plotdata = xr.concat(plotdata_test,dim='quantile')
            
            
            plt.close()
            varplot = plotdata 
            levels_plot = np.linspace(-10,10,21)
            levels_cbar = np.linspace(-10,10,11)
            plot_title  = 't2m anomaly ' + np.datetime_as_string(plotdata.validation_time[0].values, unit='D') + 'step (days)' + str(lt.days)
            fname       = 't2m_anomaly_' + np.datetime_as_string(plotdata.validation_time[0].values, unit='D') + '_step_' + str(lt.days)
            label_text  = 'K'

            im = varplot.plot(
                x               = 'lon',
                y               = 'lat',
                col              = 'quantile',
                col_wrap         = 3,
                levels           = levels_plot,
                subplot_kws      = dict(projection=ccrs.PlateCarree()),
                transform        = ccrs.PlateCarree(),
                cbar_kwargs      = {'label': label_text, 'ticks': levels_cbar}
            )


            for i,ax in enumerate(im.axes.flat):
                ax.coastlines(resolution = '10m', 
                              color      = 'black',
                              linewidth  = 0.2)
                ax.set_extent((0, 25, 55, 75), crs=ccrs.PlateCarree())
 
        
            plt.suptitle(plot_title)

            plt.savefig(fname+'.png',dpi='figure',bbox_inches='tight')
            plt.close()
            print('Figure stored at: '+fname+'.png')  
            
            
            def indices(a, func):
                return [i for (i, val) in enumerate(a) if func(val)]


            prob_dataset = []
            prob_gt   = np.empty([fcdata_sel_day.lon.shape[0],fcdata_sel_day.lat.shape[0]])
            prob_gt[:] = np.NaN
            prob_lt   = np.empty([fcdata_sel_day.lon.shape[0],fcdata_sel_day.lat.shape[0]])
            prob_lt[:] = np.NaN

            for count_lat, value_lat in enumerate(fcdata_sel_day.lat.values):
    
                for count_lon, value_lon in enumerate(fcdata_sel_day.lon.values):
    
                    temp = np.squeeze(fcdata_sel_day)
   
                    inds_gt = indices(temp[:,count_lon,count_lat], lambda x: x > 0)
                    inds_lt = indices(temp[:,count_lon,count_lat], lambda x: x < 0)
                    prob_gt[count_lon,count_lat] = (len(inds_gt)/51)*100
                    prob_lt[count_lon,count_lat] = (len(inds_lt)/51)*100
        
            prob_dataset = xr.Dataset(
                {
                    "prob_gt": (["lon", "lat"], prob_gt),
                    "prob_lt": (["lon", "lat"], prob_lt),
                },
                coords={
                    "lon": (["lon",], fcdata_sel_day.lon.values),
                    "lat": (["lat"], fcdata_sel_day.lat.values),
                },
            )
            
            plt.close()
            varplot = prob_dataset.prob_gt
            levels_plot = np.linspace(-0.1,100.1,11)
            levels_cbar = np.linspace(0,100,11)
            plot_title  = 'prob t2m anomaly > 0 ' + np.datetime_as_string(plotdata.validation_time[0].values, unit='D') + 'step (days)' + str(lt.days)
            fname       = 't2m_prob_gt_' + np.datetime_as_string(plotdata.validation_time[0].values, unit='D') + '_step_' + str(lt.days)
            label_text  = '%'

            im = varplot.plot(
                x                = 'lon',
                y                = 'lat',
                levels           = levels_plot,
                subplot_kws      = dict(projection=ccrs.PlateCarree()),
                transform        = ccrs.PlateCarree(),
                cbar_kwargs      = {'label': label_text, 'ticks': levels_cbar},
                cmap             = 'Reds',
                robust           = True
            )


            im.axes.coastlines(resolution = '10m', 
                        color      = 'black',
                        linewidth  = 0.2)
            im.axes.set_extent((0, 25, 55, 75), crs=ccrs.PlateCarree())
 
        
            plt.suptitle(plot_title)

            plt.savefig(fname+'.png',dpi='figure',bbox_inches='tight')
            plt.close()
            print('Figure stored at: '+fname+'.png')
            
            plt.close()
            varplot = prob_dataset.prob_lt
            levels_plot = np.linspace(-0.1,100.1,11)
            levels_cbar = np.linspace(0,100,11)
            plot_title  = 'prob t2m anomaly < 0 ' + np.datetime_as_string(plotdata.validation_time[0].values, unit='D') + 'step (days)' + str(lt.days)
            fname       = 't2m_prob_lt_' + np.datetime_as_string(plotdata.validation_time[0].values, unit='D') + '_step_' + str(lt.days)
            label_text  = '%'

            im = varplot.plot(
                x               = 'lon',
                y               = 'lat',
                levels           = levels_plot,
                subplot_kws      = dict(projection=ccrs.PlateCarree()),
                transform        = ccrs.PlateCarree(),
                cbar_kwargs      = {'label': label_text, 'ticks': levels_cbar},
                cmap             = 'Blues',
                robust           = True
            )


            im.axes.coastlines(resolution = '10m', 
                        color      = 'black',
                        linewidth  = 0.2)
            im.axes.set_extent((0, 25, 55, 75), crs=ccrs.PlateCarree())
 
        
            plt.suptitle(plot_title)

            plt.savefig(fname+'.png',dpi='figure',bbox_inches='tight')
            plt.close()
            print('Figure stored at: '+fname+'.png')
