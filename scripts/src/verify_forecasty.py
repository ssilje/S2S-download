import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import time as time_lib
import gridpp
import properscoring as ps
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
from S2S.local_configuration import config
from S2S.data_handler        import ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.graphics.graphics as gr
from S2S.graphics import mae,crps,brier,spread_error
import scripts.Henrik.create_domain_file

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



path_e = 't2m/'
long_name = 'absolute_t2m'
Archive().make_dir(config['VALID_DB']+path_e)

domainID = 'norwegian_coast'

var      = 't2m'

var1     = False
var2     = False

#var      = 'abs_wind'
#var1     = 'u10'
#var2     = 'v10'

t_start  = (2019,7,1)
t_end    = (2020,6,26)
# t_end    = (2020,2,3)


clim_t_start  = (1999,1,1)
clim_t_end    = (2021,1,4)

process_hindcast     = False
process_era          = False
make_time_series     = False

plot_MAE             = True

high_res             = False
verify               = True

# steps = pd.to_timedelta([7,14,23,30,37,44],'D')
steps = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')


    
###################
#### Load data ####
###################
hindcast = xr.open_dataset(config['VALID_DB']+path_e+long_name+\
                                '_hindcast.nc')[var]
hindcast_a = xr.open_dataset(config['VALID_DB']+path_e+long_name+\
                            '_anomalies_hindcast.nc')[var]
stacked_era = xr.open_dataset(config['VALID_DB']+path_e+\
                                            long_name + '_era.nc')[var]
stacked_era_a = xr.open_dataset(config['VALID_DB']+path_e+\
                                long_name + '_anomalies_era.nc')[var]
clim_mean = xr.open_dataset(config['VALID_DB']+path_e+\
                                long_name + '_era_mean.nc')[var]
clim_std = xr.open_dataset(config['VALID_DB']+path_e+\
                                long_name + '_era_std.nc')[var]
random_fc = xr.open_dataset(config['VALID_DB']+path_e+\
                                long_name + '_random-forecast.nc')[var]
random_fc_a = xr.open_dataset(config['VALID_DB']+path_e+\
                                long_name + '_random-forecast_anomalies.nc')[var]




stacked_era = xh.assign_validation_time(stacked_era)
stacked_era_a = xh.assign_validation_time(stacked_era_a)
hindcast = xh.assign_validation_time(hindcast)
hindcast_a = xh.assign_validation_time(hindcast_a)
clim_mean = xh.assign_validation_time(clim_mean)
clim_std = xh.assign_validation_time(clim_std)
random_fc = xh.assign_validation_time(random_fc)
random_fc_a = xh.assign_validation_time(random_fc_a)
    
observations = stacked_era
model = hindcast
clim_mean = clim_mean
clim_std = clim_std
dim='validation_time.month'
cc = []


for lt in steps:
#lt= steps[0]    
    mod = model.sel(step=pd.Timedelta(lt,'D'))
    obs = observations.sel(step=pd.Timedelta(lt,'D'))
    cm  = clim_mean.sel(step=pd.Timedelta(lt,'D'))
    cs  = clim_std.sel(step=pd.Timedelta(lt,'D'))
   
  
    x_group = list(mod.groupby(dim)) # lagar en liste for kvar mnd (nr_mnd, xarray)
    y_group = list(obs.groupby(dim))
    cm_group = list(cm.groupby(dim))
    cs_group = list(cs.groupby(dim))

 
    c = [] #lagar en ny xarray med score for kvar mnd
    
    for n,(xlabel,xdata) in enumerate(x_group): # loop over each validation month. n går frå 0-11, xlabel 1-12, xdata: dataene
    
        ylabel,ydata   = y_group[n]
        cmlabel,cmdata = cm_group[n]
        cslabel,csdata = cs_group[n]

        xdata  = xdata.unstack().sortby(['time']) #mod
        ydata  = ydata.unstack().sortby(['time']) # obs
        cmdata = cmdata.unstack().sortby(['time'])
        csdata = csdata.unstack().sortby(['time'])

        xdata,ydata,cmdata,csdata = xr.align(xdata,ydata,cmdata,csdata)

        score_mean   = xs.mae(xdata.mean('member',skipna=True),ydata,dim=[])
        score_clim   = xs.mae(cmdata,ydata,dim=[])
   
        SS = 1 - score_mean/score_clim
        SS = SS.median('time',skipna=True)
        #time_month = np.empty(1) 
        #print(time_month)
        
        time_month=xlabel
        #print(time_month)
        
        SS=SS.assign_coords(time_month=time_month)
        c.append(SS)
    
    skill_score = xr.concat(c,dim='time_month') ## må legge dei etter kvarandre med mnd
    skill_score_step = skill_score
    skill_score = skill_score.drop('step')
  
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
SS_step = xr.concat(cc,dim='step') 

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




