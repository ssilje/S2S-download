import pandas as pd
import xarray as xr
import xskillscore as xs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns

from S2S.data_handler import ERA5
#from S2S.process import Hindcast, Observations, Forecast, Observations_hcfc
from S2S.process import Hindcast, Observations, Forecast

from S2S.graphics import mae,crps,graphics as mae,crps,graphics
from S2S import models
import S2S.xarray_helpers    as xh
from S2S.scoring import uncentered_acc, centered_acc
from S2S.local_configuration import config


bounds          = (0,28,55,75)
var             = 'tp'

#clabel          = 'K'
t_start         = (2018,3,1)
t_end           = (2018,8,1)


clim_t_start    = (1999,1,1)
clim_t_end      = (2019,12,31)


high_res        = False

steps           = pd.to_timedelta([7, 14, 21, 28, 35, 42],'D')

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
                    
era = ERA5(high_res=high_res)\
                            .load(var,clim_t_start,clim_t_end,bounds)[var]

## need to have grid_observations_fc first. if grid_observations_hc is first, the date in grid_observations_fc is the same as grid_observations_hc.
# still something strange with hc fc obs. 

grid_observations_fc = Observations(
                                   name='Era',
                                   observations=era,
                                   forecast=grid_forecast,
                                   process=True
)
                   
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

reanalysis         = xh.assign_validation_time(grid_observations_fc.data)
forecast           = xh.assign_validation_time(grid_forecast.data)
forecast = forecast.mean(dim='member')
fcc_step = []  
re_step  = []

for lt in steps:
#lt = steps[3]
    fc_steps          = forecast.sel(step=pd.Timedelta(lt,'D')) #loop through each month
    re_steps          = reanalysis.sel(step=pd.Timedelta(lt,'D'))           

    dim               = 'validation_time.month'
    
    fc_group          = list(fc_steps.groupby(dim)) 
    re_group          = list(re_steps.groupby(dim))
    
    fcc_month = []
    re_month  = [] 



    for m,(mf,fcdata) in enumerate(fc_group): #loop through each month
        mr,eradata          = re_group[m]
    
        for reg in (
              'Norway',
            #  'MEU'
        ):
           
            fcdata_sel = fcdata.sel(lat=slice(int(region[reg]['minlat']),int(region[reg]['maxlat'])), lon=slice(int(region[reg]['minlon']),int(region[reg]['maxlon'])))
            eradata_sel = eradata.sel(lat=slice(int(region[reg]['minlat']),int(region[reg]['maxlat'])), lon=slice(int(region[reg]['minlon']),int(region[reg]['maxlon'])))
        
            fcdata_sel = fcdata_sel.mean('lon').mean('lat')
            eradata_sel = eradata_sel.mean('lon').mean('lat')
        
            fcdata_sel_df= fcdata_sel.drop('step').to_dataframe()
            eradata_sel_df = eradata_sel.drop('step').to_dataframe()
        
            plt.close()
            fig,ax=plt.subplots()
            era = ax.plot(eradata_sel_df.validation_time,eradata_sel_df.tp, color='red',label='ERA5')
            fc= ax.plot(fcdata_sel_df.validation_time,fcdata_sel_df.tp, color='blue',label='FC')
        
            ax.legend(loc='lower right')
            x_dates = eradata_sel_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
            ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
            figname = 'FC_ERA_test_' + str(lt.days) +  '_'  + str(mf) + '_' + '_month.png'
            plt.savefig(figname,dpi='figure',bbox_inches='tight')
        
        
    
