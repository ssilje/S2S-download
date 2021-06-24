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
   
   # a = mod.groupby("validation_time.month").mean().mean('member')
    x_group = list(mod.groupby(dim)) # lagar en liste for kvar mnd med 2 dim (nr_mnd, xarray)
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
        time_month = np.empty(1)
        time_month=xlabel
        
        SS=SS.assign_coords(time_month=time_month)
        c.append(SS)
    
    skill_score = xr.concat(c,dim='time_month') ## må legge dei etter kvarandre med mnd
    skill_score_step = skill_score
    skill_score = skill_score.drop('step')
    
# Bruk denne for å plotte lead time med skill    
#    im = skill_score.transpose('lat','lon','time_month').plot(
#     ...:         x='lon',
#     ...:         y='lat',
#     ...:         col='time_month',
#     ...:         col_wrap=3,
#     ...:         subplot_kws=dict(projection=ccrs.PlateCarree()),
#     ...:         transform=ccrs.PlateCarree(),
#     ...:         cmap=cmap,
#     ...:         norm=norm
#     ...:        )
    #fg,axes = plt.axes(subplot_kws=dict(projection=ccrs.PlateCarree()))
    levels = [-1,-0.5,-0.25,-0.05,0.05,0.25,0.5,1]
    im = skill_score.transpose('lat','lon','time_month').plot(
        x='lon',
        y='lat',
        col='time_month',
        col_wrap=3,
        levels=levels,
        subplot_kws=dict(projection=ccrs.PlateCarree()),
        transform=ccrs.PlateCarree(),
        cbar_kwargs={'label':'SS',
                     'ticks': levels}
    )
    
    
    for i,ax in enumerate(im.axes.flat):
        ax.coastlines(resolution='10m', color='black',\
                      linewidth=0.2)
        ax.set_title(month(i))
    

   

  
    # ax = fg.add_gridspec(3, 3)
    #cb = fg.colorbar(im, ax=[axes[-1, :]], location='bottom',boundaries=levels,extend='both') 
    #plt.tight_layout()

    #fg.suptitle('SS of MAE at lead time: '+str(lt))   
    plt.suptitle('SS of MAE at lead time: '+str(lt))
   
    plt.savefig('test_SS_day_' + str(lt.days) + '.png')

    cc.append(skill_score_step)
SS_step = xr.concat(cc,dim='step') ## må legge dei etter kvarandre med mnd
SS_group = list(SS_step.groupby('time_month'))

test = np.empty(SS_step[2,3].shape)
test[:] =np.NaN

for n,(xlabel,xdata) in enumerate(SS_group): # looping over each month
    index = xdata.where(xdata.values >0) # finding data with skill
    for lt in steps:
        skill_map= index.sel(step=pd.Timedelta(lt,'D'))
        test[:] = skill_map.where(skill_map)
        da_stacked[da_stacked.notnull()]
            
