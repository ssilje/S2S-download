import xarray as xr

import pandas as pd





def read_grib(dirbase_S2S,product,ftype,d,lat,lon):
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')
    dS2S = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,ftype,product,'.grb')
    print('reading file:')
    print(dS2S)
    dataopen = xr.open_dataset(dS2S,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    return dataopen
  
  

def read_grib_cf_pf(dirbase_S2S,product,d,lat,lon,var_short,cycle):
    dir = '%s/%s/%s/'%(dirbase_S2S,product,'/ECMWF/sfc')  
    dS2S_cf = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'cf',product,'.grb')
    dS2S_pf = '%s/%s/%s_%s_%s_%s_%s%s'%(dir,var_short,var_short,cycle,d,'pf',product,'.grb')
    
    print('reading file:')
    print(dS2S_pf)   
    dataopen_pf = xr.open_dataset(dS2S_pf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe().reset_index(level='number')
    
    print('reading file:')
    print(dS2S_cf)   
    dataopen_cf = xr.open_dataset(dS2S_cf,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() # Picking out a grid point
    
    dataopen = dataopen_cf.append(dataopen_pf).set_index('number',append=True) #merging pf and cf
    
    if product == 'forecast':
         dataopen = dataopen.set_index('time',append=True)
         dataopen.index = dataopen.index.swaplevel(1,2)
  
    return dataopen

  

def calc_stats_lead_time(dataopen,step,var,ftype):
    if ftype == 'cf':
        data_mean = dataopen.loc[(step,slice(None)),var].mean()
        data_std = dataopen.loc[(step,slice(None)),var].std()
    if ftype == 'pf':
        data_mean = dataopen.loc[(slice(None),step,slice(None)),var].mean()
        data_std = dataopen.loc[(slice(None),step,slice(None)),var].std()
        
    data_stats = pd.DataFrame({"data_mean": [data_mean], "data_std": [data_std]}, index=[step])
    data_stats.index.name = 'step'
    
    return data_stats

  
  
def calc_stats_lead_time_cf_pf(dataopen,step,var):
    data_mean = dataopen.loc[(step,slice(None)),var].mean()
    data_std = dataopen.loc[(step,slice(None)),var].std()
        
    data_stats = pd.DataFrame({"data_mean": [data_mean], "data_std": [data_std]}, index=[step])
    data_stats.index.name = 'step'
    
    return data_stats
