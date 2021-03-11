import xarray as xr
import os
import pandas as pd


def read_grib_file(
    S2S_dirbase,
    product,
    model_version,
    var_name_abbr,
    cast_type,
    date_str,
):
    file_name =  '_'.join([var_name_abbr, model_version, date_str, cast_type, product]) + '.grb'
    file_path = os.path.join(S2S_dirbase, product, 'ECMWF', 'sfc', var_name_abbr, file_name)
    
    print('reading file: ' + file_path)    
    dataopen = xr.open_dataset(file_path, engine='cfgrib')

    return dataopen



def read_grib_file_point(
    S2S_dirbase,
    product,
    model_version,
    var_name_abbr,
    cast_type,
    date_str,
    lat,
    lon,
    df_con
):
    file_name =  '_'.join([var_name_abbr, model_version, date_str, cast_type, product]) + '.grb'
    file_path = os.path.join(S2S_dirbase, product, 'ECMWF', 'sfc', var_name_abbr, file_name)
    
    print('reading file:' + file_path)
    if df_con:  
        dataopen = xr.open_dataset(file_path,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest').to_dataframe() 
    else: 
        dataopen = xr.open_dataset(file_path,engine='cfgrib').sel(latitude=lat, longitude=lon, method='nearest')
    return dataopen
  
  

def read_grib_file_merge_ftype(
    S2S_dirbase,
    product,
    model_version,
    var_name_abbr,
    date_str,
    lat,
    lon
):
    file_name_cf =  '_'.join([var_name_abbr, model_version, date_str,'cf', product]) + '.grb'
    file_path_cf = os.path.join(S2S_dirbase, product, 'ECMWF', 'sfc', var_name_abbr, file_name_cf)
    
    file_name_pf =  '_'.join([var_name_abbr, model_version, date_str,'pf', product]) + '.grb'
    file_path_pf = os.path.join(S2S_dirbase, product, 'ECMWF', 'sfc', var_name_abbr, file_name_pf)
    
   
    
    print('reading file:')
    print(file_path_pf)   
    dataopen_pf = xr.open_dataset(file_path_pf,engine='cfgrib').to_dataframe().reset_index(level='number')
    
    print('reading file:')
    print(file_path_cf)   
    dataopen_cf = xr.open_dataset(file_path_cf,engine='cfgrib').to_dataframe() # Picking out a grid point
    
    dataopen = dataopen_cf.append(dataopen_pf).set_index('number',append=True) #merging pf and cf
    
    if product == 'forecast':
         dataopen = dataopen.set_index('time',append=True)
         dataopen.index = dataopen.index.swaplevel(2,3)
  
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
