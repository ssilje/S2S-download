import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from datetime import timedelta
#import gridpp
import json 
import os




from S2S.file_handling import read_grib_file,read_grib_file_point, read_grib_file_slice_merge_ftype, check_file
from S2S.local_configuration import config

dates_monday = pd.date_range("20180503", periods=18, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20180507", periods=18, freq="7D") # forecasts start Thursday
dates_fcycle = dates_monday.union(dates_thursday)   

fc_week = {
        "week1" : ['1 days', '2 days','3 days','4 days','5 days','6 days','7 days' ],
        "week2" : ['8 days', '9 days','10 days','11 days','12 days','13 days','14 days'],
        "week3" : ['15 days', '16 days','17 days','18 days','19 days','20 days','21 days'],
        "week4" : ['22 days', '23 days','24 days','25 days','26 days','27 days','28 days']     
}
        


var_name_abbr='t2m' #tp
mdl_vrsn='CY43R3_CY45R1'

S2S_dirbase=config['S2S_DIR_summer2018']

#for d in dates_fcycle:
curr_date=dates_fcycle[0].strftime('%Y-%m-%d')
for product in (
        'forecast',
        'hindcast',
    ):            
    filecheck = check_file(
    dirbase=S2S_dirbase,
    product=product,
    model_version=mdl_vrsn,
    var_name_abbr=var_name_abbr,
    date_str=curr_date
    )
        
    if filecheck is True:
        print('file exist')
        globals()[f"ds_{product}"] = read_grib_file_slice_merge_ftype(
        dirbase=S2S_dirbase,
        product=product,
        model_version=mdl_vrsn,
        var_name_abbr=var_name_abbr,
        date_str=curr_date,
        lat=[80,45],
        lon=[0,20]
        )
    else: 
        print('Missing file')
           
        print(product)
        print(curr_date)
## de_hindcast multiindex: latitude, longitude, step, time, number    
idx_name=ds_hindcast.index.names
idx_values=ds_hindcast.index.values
# to access: idx_name[1]
ds_hindcast = ds_hindcast.drop(columns=['heightAboveGround'])
ds_forecast = ds_forecast.drop(columns=['heightAboveGround'])

LAT_unique=ds_hindcast.index.get_level_values(idx_name[0]).unique()
LON_unique=ds_hindcast.index.get_level_values(idx_name[1]).unique()

df_out_step = pd.DataFrame(
    columns = tuple([
        'data_mean', 
        'data_std',
        'step',
        'latitude', 
        'longitude' 
        
    ])
)
#df_out2.multiindex.name = ["latitude", "longitude"] #'hdate'
for ll in LAT_unique:
    for lo in LON_unique:
        for w in fc_week:
            for ww in fc_week[w]:
                data_mean=ds_hindcast.t2m[(ll,lo,ww, slice(None),slice(None))].reset_index(level='time', drop=True).mean() #mean is mean over all ens members
                data_std=ds_hindcast.t2m[(ll,lo,ww, slice(None),slice(None))].reset_index(level='time', drop=True).std() #mean is mean over all ens members
                row_dict = {
                'data_mean' : data_mean,
                'data_std' : data_std,
                'step' : ww,
                'latitude' : ll,
                'longitude' : lo,
                }
                df_out_step = df_out_step.append(pd.Series(row_dict),ignore_index=True)        
            # df_out = pd.DataFrame({"data_mean": [data_mean], "data_std": [data_std]}, index=[hdate])
           # df_out.index.name = 'hdate'

clim_hc=df_out_step.set_index('step')
clim_hc=clim_hc.set_index('latitude', append=True)
clim_hc=clim_hc.set_index('longitude', append=True)


df_out_anom = pd.DataFrame(
    columns = tuple([
        'data_anom', 
        'step',
        'latitude', 
        'longitude',
        'time',
        'number',
        'valid_time'
        
    ])
)
for ll in LAT_unique:
    for lo in LON_unique:
        for w in fc_week:
            for ww in fc_week[w]:
                 refyear = int(curr_date[:4]) 
                 for i in range(refyear-20,refyear): 
                     hdate = curr_date.replace('%i'%refyear,'%i'%i) 
                     for ens in range(1,11):  
                         data_anom=ds_hindcast.t2m[(ll,lo,ww, hdate,ens)]-clim_hc.data_mean[(ww,ll,lo)]
                         data_valid_time=ds_hindcast.valid_time[(ll,lo,ww, hdate,ens)]
                         row_dict = {
                              'data_anom' : data_anom,
                              'step' : ww,
                              'latitude' : ll,
                              'longitude' : lo,
                              'time' : hdate,
                              'number' : ens,
                              'valid_time':data_valid_time,
                              } 
                         df_out_anom = df_out_anom.append(pd.Series(row_dict),ignore_index=True)        
            
