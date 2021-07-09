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
        varplot must be a xarray.DataArray with dimension time_month (12), lon (x), lat (y)
        
    # TO DO: 
    Sjekk her: plottet blir det samme om eg brukar transpose eller ikkje. 
    #im = skill_score_at_lt.skill.transpose('lon','lat','time_month').plot(
        """
        
        
        im = varplot.plot(
                x               = 'lon',
                y               = 'lat',
               col              = 'time_month',
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
            ax.set_title(graphics.month(i))
        
        plt.suptitle(plot_title)
        plt.savefig(config['SAVEFIG']+\
                        fname+'.png',dpi='figure',bbox_inches='tight')
        plt.close()
        print('Figure stored at: '+config['SAVEFIG']+fname+'.png')
   
#graphics.save_fig(im,fname) 
#plt.savefig(plot_save)

def ACC_grid(
        forecast,
        observations,
        centered=False,
):
        
        """
        read in xarray.DataArray for: 
        - forecast with dim: time, lon, lat (use the ensemble mean or split up for each member)
        - obeservations with dim: time, lon, lat
        return 
        - ACC_dataset with dim lon, lat
        
        """
        
        
        ds = xr.merge(
                        [
                           # forecast.mean('member').rename('fc'),
                            forecast.rename('fc'),    
                            observations.rename('obs')
                        ],join='inner',compat='override'
                    )
       
      
        ACC_dataset = []
        acc    = np.empty([xdata.lon.shape[0],xdata.lat.shape[0]])
        acc[:] = np.NaN

        for nlat,ilat in enumerate(ds.lat):
            for nlon,ilon in enumerate(ds.lon): 
             
                x = ds.fc.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat],drop=True)
                y = ds.obs.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat],drop=True)
                
                if centered:
                    ACC = centered_acc(x,y)
                else:
                    ACC = uncentered_acc(x,y)

                acc[nlon,nlat] = ACC
                
        ACC_dataset = xr.Dataset(
            {
                "ACC": (["lon", "lat"], acc),
            },
            coords={
                "lon": (["lon",], xdata.lon),
                "lat": (["lat"], xdata.lat),
            },
        )
    
        return ACC_dataset
        
  
def SS_lt(
  SS_data,
):
        """
        Reads in a dataset with calculated skill-score with dim step, time_month, lat,lon
        returns a xarray.Dataset with last leadtime with skill (varname = skill) with dim lat, lon, time_month
        """
        
  
        SS_group = list(SS_data.groupby('time_month'))

        c_ss     = [] # ny xarray med siste lead time med skill 
 
        for n,(xlabel,xdata) in enumerate(SS_group): # looping over each month
    
            index       = xdata.where(xdata.values >0) # finding data with skill
    
            ss_dataset  = []
    
            test        = np.empty([xdata.lon.shape[0],xdata.lat.shape[0]])
            test[:]     = np.NaN
    
            for nlat,ilat in enumerate(xdata.lat):
                for nlon,ilon in enumerate(xdata.lon):
                
                    xdata_ss        = xdata.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat]).where(xdata.sel(lon=xdata.lon[nlon],lat=xdata.lat[nlat])>0,drop=True) # find the time steps with skill
                
                    if xdata_ss.shape[0] == 0:
                        test[nlon,nlat] = np.nan
                    else:
                        test[nlon,nlat] = xdata_ss[-1].step.dt.days.item() 
    
            ss_dataset      = xr.Dataset(
                {
                    "skill": (["lon", "lat"], test),
                },
                coords      ={
                    "lon": (["lon",], xdata.lon),
                    "lat": (["lat"], xdata.lat),
                },
            )
        
       
            ss_dataset      = ss_dataset.assign_coords(time_month=xlabel)
        
            c_ss.append(ss_dataset)
    
        skill_score_at_lt = xr.concat(c_ss,dim='time_month')  
        
        return skill_score_at_lt
  
  
  
bounds          = (0,28,55,75)

var             = 't2m'

t_start         = (2019,7,1)
t_end           = (2020,6,26)


clim_t_start    = (1999,1,1)
clim_t_end      = (2021,1,4)


high_res        = False

steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

dim             = 'validation_time.month'

cc              = []
mm              = []
Data_skill      = []
ACcc            = []

print('\tProcess Hindcast')
grid_hindcast = Hindcast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        steps=steps,
                        process=True,
                        download=False,
                        split_work=True
                    )
#    absolute:      grid_hindcast.data
#    anomalies:     grid_hindcast.data_a
#    clim_mean:     grid_hindcast.mean
#    clim_std:      grid_hindcast.std

print('\tProcess ERA5')
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



stacked_era     = xh.assign_validation_time(grid_observations.data)
stacked_era_a   = xh.assign_validation_time(grid_observations.data_a)

hindcast        = xh.assign_validation_time(grid_hindcast.data)
hindcast_a      = xh.assign_validation_time(grid_hindcast.data_a)

clim_mean = xh.assign_validation_time(grid_observations.mean)
clim_std = xh.assign_validation_time(grid_observations.std)

print('\tGenerate random forecasts')
random_fc_a = models.deterministic_gaussian_forecast(
        xr.full_like(clim_mean,0.),
        xr.full_like(clim_std,1.)
)

random_fc   = models.deterministic_gaussian_forecast(
        clim_mean,
        clim_std
)

    
observations    = stacked_era_a
model           = hindcast_a





print('\tLoop through all steps')

for lt in steps:

    mod = model.sel(step=pd.Timedelta(lt,'D'))
    obs = observations.sel(step=pd.Timedelta(lt,'D'))
    cm  = xr.full_like(observations,0).sel(step=pd.Timedelta(lt,'D')) ## er null pga bruker anomalier

    era = stacked_era.sel(step=pd.Timedelta(lt,'D'))
    hc  = hindcast.sel(step=pd.Timedelta(lt,'D'))
    
   
  
    x_group = list(mod.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    y_group = list(obs.groupby(dim))
    cm_group = list(cm.groupby(dim))

    era_group = list(era.groupby(dim))
    hc_group = list(hc.groupby(dim))
  
    m = []
    c = [] #lagar en ny xarray med score for kvar mnd
    acc = [] #lagar en ny xarray med ACC for kvar mnd
    era_tmp = []
    hc_tmp = []    

    print('\tLoop through each month')
    for n,(xlabel,xdata) in enumerate(x_group): # loop over each validation month. n går frå 0-11, xlabel 1-12, xdata: dataene
    
        ylabel,ydata   = y_group[n]
        cmlabel,cmdata = cm_group[n]
        
        eralabel,eradata = era_group[n]
        hclabel,hcdata = hc_group[n]
        
        xdata  = xdata.unstack().sortby(['time']) #mod
        ydata  = ydata.unstack().sortby(['time']) # obs
        cmdata = cmdata.unstack().sortby(['time'])
        
        eradata  = eradata.unstack().sortby(['time'])
        hcdata  = hcdata.unstack().sortby(['time'])
        

        xdata,ydata,cmdata = xr.align(xdata,ydata,cmdata)
        
        eradata,hcdata = xr.align(eradata,hcdata)
        
        
        print('\tCalculate MAE and MAESS')
        score_mean   = xs.mae(
            xdata.mean('member',skipna=True),
            ydata,
            dim=[])
        
        score_clim   = xs.mae(
            cmdata,
            ydata,
            dim=[])
   
        SS      = 1 - score_mean/score_clim
        SS      = SS.median('time',skipna=True).assign_coords(time_month=xlabel)
        c.append(SS)
        
        MAE_tmp      = score_mean.median('time',skipna=True).assign_coords(time_month=xlabel)      
        m.append(MAE_tmp)

        

        
        print('\tCalculate ACC')
        ACC_dataset = ACC_grid(
           forecast=xdata.mean('member'),
           observations=ydata,
           centered=False
        )
        
        
        ACC_dataset=ACC_dataset.assign_coords(time_month=xlabel)
        
        acc.append(ACC_dataset)
        
        
         # Calculate climatology
        era_mean = eradata.mean('time').drop('step').assign_coords(time_month=xlabel)
        era_tmp.append(era_mean)
        hc_mean  = hcdata.mean('member').mean('time').drop('step').assign_coords(time_month=xlabel)
        hc_tmp.append(hc_mean)   
        
    # Done loop month    
     
    




    skill_score = xr.concat(c,dim='time_month') ## stacking the data along month dimension
    #skill_score_step = skill_score
    #skill_score = skill_score.drop('step')
    cc.append(skill_score) # Saving MAESS for each month and step
     
    MAE = xr.concat(m,dim='time_month') ## stacking the data along month dimension
   # MAE_step_tmp = MAE 
   # MAE  = MAE.drop('step')
    mm.append(MAE) #  Saving MAE for each month and step
    
    ACC_score = xr.concat(acc,dim='time_month') 
    ACC_score = ACC_score.assign_coords(step=lt)
    ACC_score_step = ACC_score 
    ACcc.append(ACC_score_step) 

# loop leadtime done

SS_step  = xr.concat(cc,dim='step')

MAE_step  = xr.concat(mm,dim='step')

Data_skill  = xr.concat(ACcc,dim='step') 

Data_skill  = Data_skill.assign(MAESS=SS_step)

Data_skill  = Data_skill.assign(MAE=MAE_step)

Data_skill = Data_skill.assign(MAESS_best_lt=SS_lt(SS_data=SS_step).skill)

print('\tSaving calculated scores as netcdf')
outfilename = 'hindcast_skill_' + var + '.nc'
print('\t saving file with', \
      '\nMAE', Data_skill.MAE.dims,\
      '\nMAESS', Data_skill.MAESS.dims,\
      '\nMAESS_best_lt', Data_skill.MAESS_best_lt.dims,\
      '\nACC', Data_skill.MAESS_best_lt.dims)
      
Data_skill.to_netcdf(path=outfilename , mode='w')


print('\tPlotting')
## Plotting
for lt in steps:
  
    plot_months(
        varplot     = Data_skill.ACC.sel(step=lt),
        levels_plot = np.linspace(-1,1,21),
        label_text  = 'ACC',
        levels_cbar = np.linspace(-1,1,11),
        plot_title  = 'ACC',
        fname       = 'hindcast_ACC_days_' + str(lt.days) + '_' + var,
    )

    
    plot_months(
        varplot     = Data_skill.MAESS.sel(step=lt),
        levels_plot = np.linspace(-1,1,21),
        label_text  = 'MAESS',
        levels_cbar = np.linspace(-1,1,11),
        plot_title  = 'MAESS',
        fname       = 'hindcast_MAESS_days_' + str(lt.days) + '_' + var,
    )
        
    plot_months(
        varplot     = Data_skill.MAE.sel(step=lt),
        levels_plot = np.linspace(0,1,21),
        label_text  = 'MAE',
        levels_cbar = np.linspace(0,1,6),
        plot_title  = 'MAE',
        fname       = 'hindcast_MAE_days_' + str(lt.days) + '_' + var,
    )


    


plot_months(
        varplot     = Data_skill.MAESS_best_lt,
        levels_plot = [3.5, 10.5, 17.5, 24.5, 31.5, 38.5, 45.5],
        label_text  = 'lead time',
        levels_cbar = [7, 14, 21, 28, 35, 42],
        plot_title  = 'last lead time with skill',
        fname       = 'hindcast_MAESS_best_lt' + '_' + var,
    )




### PLOTTING CLIMATOLOGY 
ERA_MEAN = xr.concat(era_tmp,dim='time_month') ## stacking the data along month dimension
plot_months(
        varplot     = ERA_MEAN,
        levels_plot = np.linspace(ERA_MEAN.min(),ERA_MEAN.max(),21),
        label_text  = 'K',
        levels_cbar = np.linspace(ERA_MEAN.min(),ERA_MEAN.max(),11),
        plot_title  = 'ERA Climatology',
        fname       = 'ERA_Climatology_'  + var,
    )

HC_MEAN = xr.concat(hc_tmp,dim='time_month') ## stacking the data along month dimension
plot_months(
        varplot     = HC_MEAN,
        levels_plot = np.linspace(HC_MEAN.min(),HC_MEAN.max(),21),
        label_text  = 'K',
        levels_cbar = np.linspace(HC_MEAN.min(),HC_MEAN.max(),11),
        plot_title  = 'Hindcast Climatology',
        fname       = 'HC_Climatology_' + var,
    )
