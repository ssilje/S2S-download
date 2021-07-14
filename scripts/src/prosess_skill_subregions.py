import pandas as pd
import xarray as xr
import xskillscore as xs
import pickle

import S2S.xarray_helpers    as xh
from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

from S2S.graphics import mae,crps,graphics
from S2S import models, location_cluster

def loc(name):
    return str(location_cluster.loc_from_name(name))

var             = 't2m'
clabel          = 'K'
#var             = 'abs_wind'
t_start         = (2019,7,1)
t_end           = (2020,6,26)


clim_t_start    = (1999,1,1)
clim_t_end      = (2021,1,4)
bounds = (0,28,55,75)


high_res = False
steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')
dim             = 'validation_time.month'
grid_hindcast     = Hindcast(
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

era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]

grid_observations = Observations(
                            name='Era',
                            observations=era,
                            forecast=grid_hindcast,
                            process=False
                            )



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



# load grid for catchment
data_path = '/nird/projects/NS9853K/PROJECTS/STATKRAFT/DATA/'

with open(data_path+"catchment_boundaries_wgs84.p", "rb") as stream:
    polygon_dict = pickle.load(stream)

print('\tLoop through all steps')
CLIM_l_m_s = []
ACC_l_m_s = []
MAE_l_m_s = []

for lt in steps:

    mod         = hindcast_a.sel(step=pd.Timedelta(lt,'D'))
    obs         = stacked_era_a.sel(step=pd.Timedelta(lt,'D'))
    cm          = xr.full_like(stacked_era_a,0).sel(step=pd.Timedelta(lt,'D')) ## er null pga bruker anomalier
    obs_random  = random_fc_a.sel(step=pd.Timedelta(lt,'D'))
   

    era         = stacked_era.sel(step=pd.Timedelta(lt,'D'))
    hc          = hindcast.sel(step=pd.Timedelta(lt,'D'))
    
    x_group     = list(mod.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    y_group     = list(obs.groupby(dim))
    cm_group    = list(cm.groupby(dim))
    yr_group    = list(obs_random.groupby(dim))
  

    era_group   = list(era.groupby(dim))
    hc_group    = list(hc.groupby(dim))
    
    mae_l_m     = []
    acc_l_m     = []
    clim_l_m    =[]
    
    for n,(xlabel,xdata) in enumerate(x_group): # loop over each validation month. n går frå 0-11, xlabel 1-12, xdata: dataene
    
        ylabel,ydata     = y_group[n]
        cmlabel,cmdata   = cm_group[n]
        yrlabel,yrdata   = yr_group[n]
        
        
        eralabel,eradata = era_group[n]
        hclabel,hcdata   = hc_group[n]
        
        xdata            = xdata.unstack().sortby(['time']) #mod
        ydata            = ydata.unstack().sortby(['time']) # obs
        cmdata           = cmdata.unstack().sortby(['time'])
        yrdata           = yrdata.unstack().sortby(['time']) # random forecast
        
        eradata          = eradata.unstack().sortby(['time'])
        hcdata           = hcdata.unstack().sortby(['time'])
        

        xdata,ydata,cmdata,yrdata = xr.align(xdata,ydata,cmdata,yrdata)
        
        eradata,hcdata            = xr.align(eradata,hcdata)
        
        mae         = []
        acc         = []
        clim        = []
        
    
        for k in polygon_dict:
            print(k) 
    
            minlon,minlat,maxlon,maxlat = polygon_dict[k].bounds 
            bounds                      = (round(minlon-0.5),round(maxlon+0.5),round(minlat-0.5),round(maxlat+0.5))
    
            hindcast_sel    = xdata.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
            hindcast_sel    = hindcast_sel.assign_coords(location=k)
            hindcast_sel    = hindcast_sel.mean('lon').mean('lat')
    
            obs_sel         = ydata.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
            obs_sel         = obs_sel.assign_coords(location=k)
            obs_sel         = obs_sel.mean('lon').mean('lat')
        
            clim_sel        = cmdata.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
            clim_sel        = clim_sel.assign_coords(location=k)
            clim_sel        = clim_sel.mean('lon').mean('lat')
        
            rf_sel          = yrdata.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
            rf_sel          = rf_sel.assign_coords(location=k)
            rf_sel          = rf_sel.mean('lon').mean('lat')
        
            eraclim_sel    = eradata.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
            eraclim_sel    = eraclim_sel.assign_coords(location=k)
            eraclim_sel    = eraclim_sel.mean('lon').mean('lat')
    
            hcclim_sel         = hcdata.sel(lon=slice(bounds[0],bounds[1]),lat=slice(bounds[2],bounds[3]))
            hcclim_sel         = hcclim_sel.assign_coords(location=k)
            hcclim_sel         = hcclim_sel.mean('lon').mean('lat')
    
    
            score_mean   = xs.mae(
                hindcast_sel.mean('member',skipna=True),
                obs_sel,
                dim=[])
    
            score_clim   = xs.mae(
                clim_sel,
                obs_sel,
                dim=[])
   
            SS_mae      = 1 - score_mean/score_clim
    
            MAE_dataset = xr.merge(
                                 [
                                     score_mean.median('time',skipna=True).rename('MAE'),
                                     SS_mae.median('time',skipna=True).rename('MAESS')
                                 ],join='inner',compat='override'
                             )
            MAE_dataset      = MAE_dataset.assign_coords(time_month=xlabel)      
            mae.append(MAE_dataset)
    
    
            ACC_hc = uncentered_acc(hindcast_sel.mean('member'),obs_sel) 
            ACC_rf = uncentered_acc(rf_sel,obs_sel) 
        
            ACC_dataset = xr.merge(
                                 [
                                  ACC_hc.rename('ACC_hc'),
                                  ACC_rf.rename('ACC_rf')
                               ],join='inner',compat='override'
                             )
            ACC_dataset     = ACC_dataset.assign_coords(time_month=xlabel)      
            acc.append(ACC_dataset)
        
            eraclim_sel_mean = eraclim_sel.mean('time').assign_coords(time_month=xlabel) 
            hcclim_sel_mean  = hcclim_sel.mean('member').mean('time').assign_coords(time_month=xlabel)
            
            
            CLIM_dataset = xr.merge(
                                 [
                                  eraclim_sel_mean.rename('ClimERA'),
                                  hcclim_sel_mean.rename('ClimHC')
                               ],join='inner',compat='override'
                             )
            
            clim.append(CLIM_dataset)
            
        # Done loop catchment   
     
        mae_l   = xr.concat(mae,dim='location') 
        mae_l_m.append(mae_l)
        acc_l   = xr.concat(acc,dim='location') 
        acc_l_m.append(acc_l)
        clim_l  = xr.concat(clim,dim='location') 
        clim_l_m.append(clim_l)
    # Done loop month 
    mae_l_m_s   = xr.concat(mae_l_m,dim='time_month') 
    MAE_l_m_s.append(mae_l_m_s)
    
    acc_l_m_s   = xr.concat(acc_l_m,dim='time_month') 
    ACC_l_m_s.append(acc_l_m_s)
    
    clim_l_m_s   = xr.concat(clim_l_m,dim='time_month') 
    CLIM_l_m_s.append(clim_l_m_s)
 # Done loop step 
MAE   = xr.concat(MAE_l_m_s,dim='step')
ACC   = xr.concat(ACC_l_m_s,dim='step')
CLIM  = xr.concat(CLIM_l_m_s,dim='step')
            
