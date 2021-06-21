import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import time as time_lib
import gridpp
import properscoring as ps
from S2S.local_configuration import config
from S2S.data_handler        import ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.graphics.graphics as gr

import scripts.Henrik.create_domain_file

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
    
    
if verify:
    val_obs = stacked_era
    clim_fc  = models.clim_fc(clim_mean,clim_std)
    hindcast = hindcast.sel(time=slice(val_obs.time.min(),val_obs.time.max()))
    
    

    crps_mod  = ps.crps_ensemble(np.moveaxis(val_obs.values,0,1),np.moveaxis(hindcast.values,0,-1))
    
    #crps_mod  = xs.crps_ensemble(np.moveaxis(val_obs.values,0,1),hindcast.values,dim=[]) fekke ein feilmelding
    crps_clim      = ps.crps_gaussian(x=sst_obs, mu=0, sig=std_obs) # does not weight with cos(lat), might be problematic in next operation crps_SS
    crps_clim = ps.crps_gaussian(
                                np.moveaxis(stacked_era_a.values,0,1),
                                np.moveaxis(clim_mean.values,-1,1),
                                np.moveaxis(clim_std.values,-1,1),
                                dim=[]
                                )

    gr.skill_plot(crps_mod,crps_clim,title='rawEC',filename='rawEC')
    gr.qq_plot(val_obs[var],hindcast[var],y_axlabel='raw EC',filename='rawEC')


stacked_era = xh.assign_validation_time(
                    stacked_era.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

stacked_era_a = xh.assign_validation_time(
                    stacked_era_a.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

hindcast = xh.assign_validation_time(
                    hindcast.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

hindcast_a = xh.assign_validation_time(
                    hindcast_a.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

clim_mean = xh.assign_validation_time(
                    clim_mean.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

clim_std = xh.assign_validation_time(
                    clim_std.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

random_fc = xh.assign_validation_time(
                    random_fc.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

random_fc_a = xh.assign_validation_time(
                    random_fc_a.isel(lon=5,lat=5).expand_dims('location')\
                        .assign_coords(location=['lon: 63. lat: 7.5'])
                        )

    
gr.timeseries(
                    stacked_era,
                    cast=[random_fc],
                    title='EC',
                    filename=long_name + '_random',
                    clabs=['random_fc'],
                    lead_time=[7, 14]
                )
