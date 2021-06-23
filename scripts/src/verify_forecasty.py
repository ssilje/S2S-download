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

for lt in steps:
#lt= steps[0]    
    mod = model.sel(step=pd.Timedelta(lt,'D'))
    obs = observations.sel(step=pd.Timedelta(lt,'D'))
    cm  = clim_mean.sel(step=pd.Timedelta(lt,'D'))
    cs  = clim_std.sel(step=pd.Timedelta(lt,'D'))
   
    dim='validation_time.month'
    x_group = list(mod.groupby(dim))
    y_group = list(obs.groupby(dim))
    cm_group = list(cm.groupby(dim))
    cs_group = list(cs.groupby(dim))

    fg,axes = plt.subplots(ncols=4,nrows=3,\
                    subplot_kw=dict(projection=ccrs.PlateCarree()))
    axes_f = axes.flatten()
    cmap   = mpl.colors.ListedColormap(
                    ['red','red','red','white','lightblue','royalblue','blue']
                                        )
    levels = [-1,-0.5,-0.25,-0.05,0.05,0.25,0.5,1]
    norm   = BoundaryNorm(levels,cmap.N)
    for n,(xlabel,xdata) in enumerate(x_group): # loop over each validation month
    
        ylabel,ydata   = y_group[n]
        cmlabel,cmdata = cm_group[n]
        cslabel,csdata = cs_group[n]

        xdata  = mod.unstack().sortby(['time']) #mod
        ydata  = obs.unstack().sortby(['time']) # obs
        cmdata = cm.unstack().sortby(['time'])
        csdata = cs.unstack().sortby(['time'])

        xdata,ydata,cmdata,csdata = xr.align(xdata,ydata,cmdata,csdata)

        score_mean   = xs.mae(xdata.mean('member',skipna=True),ydata,dim=[])
        score_clim   = xs.mae(cmdata,ydata,dim=[])
   
        SS = 1 - score_mean/score_clim
        SS = SS.median('time',skipna=True)
    
    

        im = SS.transpose('lat','lon').plot(
            ax=axes_f[n],
            transform=ccrs.PlateCarree(),  # this is important!
            add_colorbar=False, 
            cmap=cmap, 
            norm=norm
        )
    
        ax = axes_f[n]
        

        ax.coastlines(resolution='10m', color='black',\
                     linewidth=0.2)
                                
    
        ax.set_title(month(xlabel))
        ax = fg.add_gridspec(3, 3)
    cb = fg.colorbar(im, ax=[axes[-1, :]], location='bottom',boundaries=levels) 
    plt.tight_layout()

    fg.suptitle('SS of MAE at lead time: '+str(lt))    
    fg.savefig('test_SS_day_' + str(lt.days) + '.png')

