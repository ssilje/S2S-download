import pandas as pd
import xarray as xr
import xskillscore as xs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from S2S.data_handler import ERA5
from S2S.process import Hindcast, Observations

from S2S.graphics import mae,crps,graphics as mae,crps,graphics

from S2S import models

import S2S.xarray_helpers    as xh

from S2S.scoring import uncentered_acc, centered_acc

from S2S.local_configuration import config

## Routines that are used which can be moved to S2S-folder

def plot_months(
        varplot,
        levels_plot,
        label_text,
        levels_cbar,
        plot_title,
        fname,
):
        """
        # # PLOTTING # # 
        Generates a figure with 12 subfigures (each month)
        
    # Sjekk her: plottet blir det samme om eg brukar transpose eller ikkje. 
    #im = skill_score_at_lt.skill.transpose('lon','lat','time_month').plot(
        """
        
        
    im = varplot.plot(
            x               = 'lon',
            y               = 'lat',
            col             = 'time_month',
            col_wrap        = 3,
            levels          = levels_plot,
            subplot_kws     = dict(projection=ccrs.PlateCarree()),
            transform       = ccrs.PlateCarree(),
            cbar_kwargs     = {'label': label_text, 'ticks': levels_cbar}
    )
  
    for i,ax in enumerate(im.axes.flat):
        ax.coastlines(resolution = '10m', 
                      color      = 'black',
                      linewidth  = 0.2)
        ax.set_title(graphics.month(i))
        
    plt.suptitle(plot_title)
    plt.savefig(config['SAVEFIG']+\
                    fname+'.png',dpi='figure',bbox_inches='tight')
    plt.close()
    print('Figure stored at: '+config['SAVEFIG']+filename+'.png')
   
#graphics.save_fig(im,fname) 
#plt.savefig(plot_save)

def ACC_grid(
        forecast,
        observations,
):
    
    ds = xr.merge(
                    [
                        forecast.mean('member').rename('fc'),
                        observations.rename('obs')
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
    
    return ACC_dataset
        
  
def SS_best_lt(
  SS_data,
):
   # SS_data must be a dataset with calculated skill-score for lat,lon for each month and step
  SS_group = list(SS_data.groupby('time_month'))

  c_ss =[] # ny xarray med siste lead time med skill 
 
  for n,(xlabel,xdata) in enumerate(SS_group): # looping over each month
    
      index = xdata.where(xdata.values >0) # finding data with skill
    
      ss_dataset = []
    
      test = np.empty([xdata.lon.shape[0],xdata.lat.shape[0]])
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
  return skill_score_at_lt
  
  
  
bounds = (0,28,55,75)
var      = 't2m'

t_start  = (2019,7,1)
t_end    = (2020,6,26)
#t_end    = (2019,8,29)

clim_t_start  = (1999,1,1)
clim_t_end    = (2021,1,4)

plot_MAE = True
high_res = False
steps = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

dim='validation_time.month'
cc = []
Data_skill = []

ACcc = []


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



stacked_era = xh.assign_validation_time(grid_observations.data)
stacked_era_a = xh.assign_validation_time(grid_observations.data_a)
hindcast = xh.assign_validation_time(grid_hindcast.data)
hindcast_a = xh.assign_validation_time(grid_hindcast.data_a)
clim_mean = xh.assign_validation_time(grid_observations.mean)
clim_std = xh.assign_validation_time(grid_observations.std)

    
observations = stacked_era_a
model = hindcast_a
clim_mean = xr.full_like(observations,0)



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
        
        
        SS=SS.assign_coords(time_month=time_month)
        c.append(SS)

        # Calculate ACC
        ACC_dataset = ACC_grid(
           forecast=xdata,
           observations=ydata,
        )
  
        ACC_dataset=ACC_dataset.assign_coords(time_month=time_month)
        
        acc.append(ACC_dataset)
        
    skill_score = xr.concat(c,dim='time_month') ## må legge dei etter kvarandre med mnd
    skill_score_step = skill_score
    skill_score = skill_score.drop('step')
    cc.append(skill_score_step) # lagrar MAE for kvar mnd og kvar step
    
    ACC_score = xr.concat(acc,dim='time_month') 
    ACC_score = ACC_score.assign_coords(step=lt)
    ACC_score_step = ACC_score 
    ACcc.append(ACC_score_step) 

# loop leadtime done

SS_step  = xr.concat(cc,dim='step')
Data_skill  = xr.concat(ACcc,dim='step') 
Data_skill  = Data_skill.assign(MAESS=SS_step)

Data_skill = Data_skill.assign(MAESS_best_lt=SS_best_lt(SS_data=SS_step).skill)

outfilename = 'hindcast_skill_' + var + '.nc'
Data_skill.to_netcdf(path=outfilename , mode='w')



## Plotting
for lt in steps:
  
    plot_months(
        varplot     = Data_skill.ACC.sel(step=lt),
        levels_plot = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1],
        label_text  = 'ACC',
        levels_cbar = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1],
        plot_title  = 'ACC',
        fname       = 'hindcast_ACC_days_' + str(lt.days),
    )

    
    plot_months(
        varplot     = Data_skill.MAESS.sel(step=lt),
        levels_plot = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1],
        label_text  = 'MAESS',
        levels_cbar = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1],
        plot_title  = 'MAESS',
        fname       = 'hindcast_MAESS_days_' + str(lt.days),
    )
    
plot_months(
        varplot     = Data_skill.MAESS_best_lt,
        levels_plot = [3.5, 10.5, 17.5, 24.5, 31.5, 38.5, 45.5],
        label_text  = 'lead time',
        levels_cbar = [7, 14, 21, 28, 35, 42],
        plot_title  = 'last lead time with skill',
        fname       = 'hindcast_MAESS_best_lt',
    )
