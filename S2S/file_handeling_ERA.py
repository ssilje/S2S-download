import xarray as xr
import pandas as pd

def read_ERA5_timeseries(
    dirbase,
    var_long,
    start_date,
    end_date,
    lat,
    lon,
    daymean,
):
    if end_date is False:
        
        file = '%s/%s_%s%s'%(dirbase,var_long,start_date,'.nc')
        dataopen = xr.open_dataset(file)
        if daymean is True:
            ERA5 = dataopen.sel(latitude=slice(lat[0],lat[1]),longitude=slice(lon[0],lon[1])).resample(time='D').mean()
        else:
            ERA5 = dataopen.sel(latitude=slice(lat[0],lat[1]),longitude=slice(lon[0],lon[1]))
    else:
        date=pd.date_range(start_date, end_date)
        for d in date:
            file = '%s/%s_%s%s'%(dirbase,var_long,d.strftime('%Y%m%d'),'.nc')
            dataopen = xr.open_dataset(file)
            if daymean is True:
                data_df = dataopen.sel(latitude=slice(lat[0],lat[1]),longitude=slice(lon[0],lon[1])).resample(time='D').mean().to_dataframe()
            else:
                data_df = dataopen.sel(latitude=slice(lat[0],lat[1]),longitude=slice(lon[0],lon[1])).to_dataframe()
        
            if d == date[0]:
                ERA5_df = data_df
            else: 
                ERA5_df=ERA5_df.append(data_df)
            ERA5=ERA5_df.to_xarray() 
        
    return ERA5

def read_ERA5_clim_anom(
    dirbase,
    var_long,
    date,
    lat,
    lon,
):
    for i in range(1,21):
        curr_date=dates_fcycle[0] - pd.DateOffset(years=i)
        file = '%s/%s_%s%s'%(dirbase,var_long,curr_date.strftime('%Y%m%d'),'.nc')
        dataopen = xr.open_dataset(file)
        data_df = dataopen.sel(latitude=slice(lat[0],lat[1]),longitude=slice(lon[0],lon[1])).resample(time='D').mean().to_dataframe()
        if i == 1:
            ERA5_hc_df = data_df
        else: 
            ERA5_hc_df=ERA5_hc_df.append(data_df)
        ERA5_hc=ERA5_hc_df.to_xarray() 
        
        return ERA5
        
    
