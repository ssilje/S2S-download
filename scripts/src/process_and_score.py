import pandas as pd
import xarray as xr
import xskillscore as xs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models

import S2S.xarray_helpers    as xh

from S2S.scoring import uncentered_acc, centered_acc, ACc





def month(ii):
    """
    ii : integer
    returns
    month string_tools
    """
    try:
        ii = int(ii)
        ii = ['JAN','FEB','MAR',\
            'APR','MAI','JUN',\
            'JUL','AUG','SEP',\
            'OCT','NOV','DES'][ii-1]
    except (ValueError,IndexError) as error:
        pass
    return ii


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
        subplot_kws=dict(projection=ccrs.PlateCarree()),
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
    
bounds = (0,28,55,75)
#var      = 't2m'
var      = 'abs_wind'

t_start  = (2019,7,1)
t_end    = (2020,6,26)
#t_end    = (2019,8,29)

clim_t_start  = (1999,1,1)
clim_t_end    = (2021,1,4)

plot_MAE = True
high_res = False
steps = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

grid_hindcast = Hindcast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        steps=steps,
                        process=False,
                        download=False,
                        split_work=True
                    )
#    absolute:      grid_hindcast.data
#    anomalies:     grid_hindcast.data_a
#    clim_mean:     grid_hindcast.mean
#    clim_std:      grid_hindcast.std



print('\tLoad ERA')
if var == 'abs_wind':
    era_u = ERA5(high_res=high_res)\
                         .load('u10',clim_t_start,clim_t_end,domainID)['u10]

    era_v = ERA5(high_res=high_res)\
                        .load('v10',clim_t_start,clim_t_end,domainID)['v10]

    era_u,era_v = xr.align(era_u,era_v)

    era = xr.apply_ufunc(
                    xh.absolute,era_u,era_v,
                    input_core_dims  = [[],[]],
                    output_core_dims = [[]],
                    vectorize=True,dask='parallelized'
                )

    era = era.rename(var)

else:

    era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]

grid_observations = Observations(
                            name='Era',
                            observations=era,
                            forecast=grid_hindcast,
                            process=False
                            )


#    absolute:      grid_observationst.data
#    anomalies:     grid_observationst.data_a
#    clim_mean:     grid_observationst.mean
#    clim_std:      grid_observationst.std

print('\tGenerate random forecasts')
random_fc_a = models.deterministic_gaussian_forecast(
                                            xr.full_like(grid_observations.mean,0.),
                                            xr.full_like(grid_observations.std,1.)
                                            )
random_fc   = models.deterministic_gaussian_forecast(
                                            grid_observations.mean,
                                            grid_observations.std
                                            )


stacked_era = xh.assign_validation_time(grid_observations.data)
stacked_era_a = xh.assign_validation_time(grid_observations.data_a)
hindcast = xh.assign_validation_time(grid_hindcast.data)
hindcast_a = xh.assign_validation_time(grid_hindcast.data_a)
clim_mean = xh.assign_validation_time(grid_observations.mean)
clim_std = xh.assign_validation_time(grid_observations.std)
random_fc = xh.assign_validation_time(random_fc)
random_fc_a = xh.assign_validation_time(random_fc_a)
    
observations = stacked_era_a
model = hindcast_a
clim_mean = xr.full_like(observations,0)



dim='validation_time.month'
cc = []
ACcc = []





for lt in steps:

    mod = model.sel(step=pd.Timedelta(lt,'D'))
    obs = observations.sel(step=pd.Timedelta(lt,'D'))
    cm  = clim_mean.sel(step=pd.Timedelta(lt,'D'))
    
   
  
    x_group = list(mod.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    y_group = list(obs.groupby(dim))
    cm_group = list(cm.groupby(dim))
  

 
    c = [] #lagar en ny xarray med score for kvar mnd
    acc = [] #lagar en ny xarray med ACC for kvar mnd
    
    for n,(xlabel,xdata) in enumerate(x_group): # loop over each validation month. n går frå 0-11, xlabel 1-12, xdata: dataene
    
        ylabel,ydata   = y_group[n]
        cmlabel,cmdata = cm_group[n]
        
        xdata  = xdata.unstack().sortby(['time']) #mod
        ydata  = ydata.unstack().sortby(['time']) # obs
        cmdata = cmdata.unstack().sortby(['time'])
        

        xdata,ydata,cmdata = xr.align(xdata,ydata,cmdata)
       
    # Calculate MAE
        score_mean   = xs.mae(
            xdata.mean('member',skipna=True),
            ydata,
            dim=[])
        score_clim   = xs.mae(
            cmdata,
            ydata,
            dim=[])
   
        SS = 1 - score_mean/score_clim
    
        SS = SS.median('time',skipna=True)
        time_month=xlabel
        #print(time_month)
        
        SS=SS.assign_coords(time_month=time_month)
        c.append(SS)
        
      # Calculate ACC
       
        ds = xr.merge(
                    [
                        xdata.mean('member').rename('fc'),
                        ydata.rename('obs')
                    ],join='inner',compat='override'
                )
        ACC_dataset = []
        acc_test = np.empty([xdata.lon.shape[0],xdata.lat.shape[0]])
        acc_test[:] =np.NaN
        for nlat,ilat in enumerate(ds.lat):
            for nlon,ilon in enumerate(ds.lon): 
                x = ds.fc.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat],drop=True)
                y = ds.obs.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat],drop=True)

               
                ACC = uncentered_acc(x,y)
                acc_test[nlon,nlat] = ACC
                
        ACC_dataset = xr.Dataset(
            {
                "ACC": (["lon", "lat"], acc_test),
            },
            coords={
                "lon": (["lon",], xdata.lon),
                "lat": (["lat"], xdata.lat),
            },
        )
        time_month=xlabel
        ACC_dataset=ACC_dataset.assign_coords(time_month=time_month)
        
        acc.append(ACC_dataset)
    
    
    
    skill_score = xr.concat(c,dim='time_month') ## må legge dei etter kvarandre med mnd
    skill_score_step = skill_score
    skill_score = skill_score.drop('step')
    
    ACC_score = xr.concat(acc,dim='time_month') 
    ACC_score = ACC_score.assign_coords(step=lt)
    ACC_score_step = ACC_score 
    if plot_MAE:
        
        plot_months(
            varplot = skill_score.transpose('lat','lon','time_month'),
            levels_plot = [-1,-0.5,-0.25,-0.05,0.05,0.25,0.5,1],
            label_text  = 'SS',
            levels_cbar = [-1,-0.5,-0.25,-0.05,0.05,0.25,0.5,1],
            plot_title  = 'SS of MAE at lead time: '+str(lt),
            plot_save   = 'test_SS_day_' + str(lt.days) + '.png',
        )
 

    cc.append(skill_score_step) # lagrar MAE for kvar mnd og kvar step
    ACcc.append(ACC_score_step) 
    
SS_step = xr.concat(cc,dim='step')
ACC_step = xr.concat(ACcc,dim='step') 
ACC_step.assign(MAESS=SS_step)


SS_group = list(SS_step.groupby('time_month'))

c_ss =[] # ny xarray med siste lead time med skill 

for n,(xlabel,xdata) in enumerate(SS_group): # looping over each month
    
    index = xdata.where(xdata.values >0) # finding data with skill
    
    ss_dataset = []
    
    test = np.empty(SS_step[2,3].shape)
    test[:] =np.NaN
    
    for nlat,ilat in enumerate(xdata.lat):
        for nlon,ilon in enumerate(xdata.lon):
            xdata_ss = xdata.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat]).where(xdata.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat])>0,drop=True) # find the time steps with skill
            if xdata_ss.shape[0] == 0:
                test[nlon,nlat] = np.nan
            else:
                test[nlon,nlat] = xdata_ss[-1].step.dt.days.item() 
    
    ss_dataset = xr.Dataset(
        {
            "skill": (["lon", "lat"], test),
        },
        coords={
            "lon": (["lon",], xdata.lon),
            "lat": (["lat"], xdata.lat),
        },
    )
    
    time_month=xlabel
    ss_dataset=ss_dataset.assign_coords(time_month=time_month)
    c_ss.append(ss_dataset)
    
skill_score_at_lt = xr.concat(c_ss,dim='time_month')     

plot_months(
    varplot = skill_score_at_lt.skill,
    levels_plot = [3.5, 10.5, 17.5, 24.5, 31.5, 38.5, 45.5],
    label_text  = 'lead time',
    levels_cbar = [7, 14, 21, 28, 35, 42],
    plot_title  = 'last lead time with skill',
    plot_save   = 'test_SS_leadtime.png',
)
