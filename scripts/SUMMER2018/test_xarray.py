#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
#import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from datetime import timedelta

import json 
import os




from S2S.file_handling import read_grib_file,read_grib_file_point, read_grib_file_slice_merge_ftype, check_file,read_grib_slice_mft_xarray
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
        globals()[f"xs_{product}"] = read_grib_slice_mft_xarray(
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
        
# t2m = xs_hindcast.sel(step=slice('1 days', '2 days')).t2m
#t2m = xs_hindcast.sel(step=slice('1 days', '2 days')).mean(dim='step').t2m


## reading ERA
#xs_hindcast.valid_time.sel(number=0).to_dataframe().sort_values(by='valid_time') Kan kanskje bruke denne til å få dei forskjellige datoane frå s2s

        
ERA5 = read_ERA5_timeseries(
        dirbase=config['ERA5_daily_DIR'],
        var_long='2m_temperature',
        start_date='19980501',
        end_date='20181001',
        lat=[80,45],
        lon=[0,20],
        daymean=True,
)






clim_mean = xs_hindcast.mean(dim='number').mean(dim='time')      # does this give mean over all years?
clim_std = xs_hindcast.std(dim='number').std(dim='time')      # does this give mean over all years?
anomaly_hc = xs_hindcast- clim_mean
anomaly_fc = xs_forecast- clim_mean
anom_reg_hc = anomaly_hc.mean(dim='latitude').mean(dim='longitude')
anom_reg_fc = anomaly_fc.mean(dim='latitude').mean(dim='longitude')


anom_week1_hc=anom_reg_hc.t2m.sel(step=slice('1 days', '7 days')).mean(dim='step').to_dataframe()
anom_week2_hc=anom_reg_hc.t2m.sel(step=slice('8 days', '14 days')).mean(dim='step').to_dataframe()
anom_week3_hc=anom_reg_hc.t2m.sel(step=slice('15 days', '21 days')).mean(dim='step').to_dataframe()
anom_week4_hc=anom_reg_hc.t2m.sel(step=slice('22 days', '28 days')).mean(dim='step').to_dataframe()

anom_week1_fc=anom_reg_fc.t2m.sel(step=slice('1 days', '7 days')).mean(dim='step').to_dataframe()
anom_week2_fc=anom_reg_fc.t2m.sel(step=slice('8 days', '14 days')).mean(dim='step').to_dataframe()
anom_week3_fc=anom_reg_fc.t2m.sel(step=slice('15 days', '21 days')).mean(dim='step').to_dataframe()
anom_week4_fc=anom_reg_fc.t2m.sel(step=slice('22 days', '28 days')).mean(dim='step').to_dataframe()



data_week1_hc = pd.DataFrame(anom_week1_hc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week1-hc"})
data_week2_hc = pd.DataFrame(anom_week2_hc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week2-hc"})
data_week3_hc = pd.DataFrame(anom_week3_hc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week3-hc"})
data_week4_hc = pd.DataFrame(anom_week4_hc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week4-hc"})

data_week1_fc = pd.DataFrame(anom_week1_fc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week1-fc"})
data_week2_fc = pd.DataFrame(anom_week2_fc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week2-fc"})
data_week3_fc = pd.DataFrame(anom_week3_fc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week3-fc"})
data_week4_fc = pd.DataFrame(anom_week4_fc.t2m.reset_index(level='time', drop=True).reset_index(level='number', drop=True)).rename(columns={"t2m":"t2m-week4-fc"})

week_data = pd.concat([data_week1_hc, data_week1_fc, data_week2_hc,data_week2_fc,data_week3_hc,data_week3_fc,data_week4_hc,data_week4_fc])
#fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
#boxplot = week_data.boxplot(column=['t2m-week1-hc',
 #                                   't2m-week1-fc',
 #                                   't2m-week2-hc',
 ##                                   't2m-week2-fc',
 #                                   't2m-week3-hc',
 #                                   't2m-week3-fc',
 #                                   't2m-week4-hc',
 #                                   't2m-week4-fc']
 #                           ,rot=90, fontsize=15)



#fig.savefig('test.png')    
 
        
        
        #data_anome_weeks = data_anome_weeks.append(pd.Series(row_dict), ignore_index = True)
#data_anome_weeks = data_anome_weeks.concat([pd.Series(row_dict)], ignore_index = True, axis=1)
#data_anome_weeks = data_anome_weeks.append(row_dict, ignore_index = True)
        
#boxplot = anom_week1_hc.boxplot(column=['t2m'])
#boxplot = data_anome_weeks.boxplot(column=['t2m-week1-hc',
#                                         't2m-week1-fc',
#                                         't2m-week2-hc',
#                                         't2m-week2-fc',
#                                         't2m-week3-hc',
#                                         't2m-week3-fc',
#                                         't2m-week4-hc',
#                                         't2m-week4-fc']) 



#data_anom_hc=xs_hindcast.sel(step=ww,time=hdate,number=ens)-clim_mean.sel(step=ww)
#data_anom_hc.sel(latitude=slice(65,60),longitude=slice(5,10))
                
                                    
                       
                         
   
