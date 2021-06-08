#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from datetime import timedelta
import seaborn as sns

import json 
import os




from S2S.file_handling import read_grib_file,read_grib_file_point, read_grib_file_slice_merge_ftype, check_file,read_grib_slice_mft_xarray
from S2S.local_configuration import config
from S2S.file_handeling_ERA import read_ERA5_clim_anom

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

S2S_dirbase='/nird/projects/NS9853K/DATA/S2S/SUMMER2018/' 
#config['S2S_DIR_summer2018']
print(S2S_dirbase)
for d in dates_fcycle:
    #curr_date=dates_fcycle[0].strftime('%Y-%m-%d')
    curr_date=d.strftime('%Y-%m-%d')
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

    #dates = pd.date_range("20180503", periods=28) # forecasts start Thursday
    dates = pd.date_range(curr_date, periods=28) # forecasts start Thursday
    for d in dates:
            ERA5 = read_ERA5_clim_anom(
                #dirbase=config['ERA5_daily_DIR'],
                    dirbase='/nird/projects/NS9853K/DATA/SFE/ERA_daily_nc/',
                    var_long='2m_temperature',
                    date=d,
                    lat=[80,45],
                    lon=[0,20],
            )
            ERA5_anom = ERA5.mean(dim='latitude').mean(dim='longitude').to_dataframe()
            print(ERA5_anom)
            if d == dates[0]:
                ERA5_anom_df = ERA5_anom
            else: 
                ERA5_anom_df= ERA5_anom_df.append(ERA5_anom)
    
    date_step_week1=xs_forecast.sel(
            step=slice('1 days', '7 days'),
            number=0
    ).valid_time.mean(
            dim='latitude'
    ).mean(
            dim='longitude'
    ).to_dataframe().drop(
            columns=['number']
    ).reset_index(
            level=1, drop=True
    ) 


    date_step_week2=xs_forecast.sel(
            step=slice('8 days', '14 days'),
            number=0
    ).valid_time.mean(
            dim='latitude'
    ).mean(
            dim='longitude'
    ).to_dataframe().drop(
            columns=['number']
    ).reset_index(
            level=1, drop=True
    ) 

    date_step_week3=xs_forecast.sel(
            step=slice('15 days', '21 days'),
            number=0
    ).valid_time.mean(
            dim='latitude'
    ).mean(
            dim='longitude'
    ).to_dataframe().drop(
            columns=['number']
    ).reset_index(
            level=1, drop=True
    ) 

    date_step_week4=xs_forecast.sel(
            step=slice('22 days', '28 days'),
            number=0
    ).valid_time.mean(
            dim='latitude'
    ).mean(
            dim='longitude'
    ).to_dataframe().drop(
            columns=['number']
    ).reset_index(
            level=1, drop=True
    ) 

    mean_ERA5_week1=ERA5_anom_df.loc[date_step_week1.valid_time[0].strftime('%Y%m%d'):date_step_week1.valid_time[-1].strftime('%Y%m%d')].mean()
    mean_ERA5_week2=ERA5_anom_df.loc[date_step_week2.valid_time[0].strftime('%Y%m%d'):date_step_week2.valid_time[-1].strftime('%Y%m%d')].mean()
    mean_ERA5_week3=ERA5_anom_df.loc[date_step_week3.valid_time[0].strftime('%Y%m%d'):date_step_week3.valid_time[-1].strftime('%Y%m%d')].mean()
    mean_ERA5_week4=ERA5_anom_df.loc[date_step_week4.valid_time[0].strftime('%Y%m%d'):date_step_week4.valid_time[-1].strftime('%Y%m%d')].mean()

    mean_ERA5_df = pd.DataFrame({"t2m-mean":[mean_ERA5_week1.t2m,mean_ERA5_week2.t2m,mean_ERA5_week3.t2m,mean_ERA5_week4.t2m]},index=['week1','week2','week3','week4'])
    mean_ERA5_df = pd.DataFrame({
            "t2m-week1-hc":np.nan,
            "t2m-week1-fc":mean_ERA5_week1.t2m,
            "t2m-week12-hc":np.nan,
            "t2m-week2-fc":mean_ERA5_week2.t2m,
            "t2m-week3-hc":np.nan,
            't2m-week3-fc':mean_ERA5_week3.t2m,
            "t2m-week4-hc":np.nan,
            't2m-week4-fc':mean_ERA5_week4.t2m
    }, index = ['t2m'])
    mean_ERA5_df.to_csv('era5_anom.csv')


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

#era_anom=mean_ERA5_df.iloc[:,1:9]

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
    bplot=sns.boxplot(data=week_data, 
                     width=0.5,
                     palette="colorblind",ax=axes).set(title='forecast initialized ' + curr_date)
 
    bplot=sns.stripplot(data=mean_ERA5_df, 
                       jitter=True, 
                       size=20, 
                       marker="D",
                       alpha=0.5,
                       color='black',ax=axes)
    bplot.set_xticklabels(['HC w1','FC w1','HC w2','FC w2','HC w3','FC w3','HC w4','HC w4'],rotation=90)

    fig.savefig(curr_date + '_t2manom.png')    
 

