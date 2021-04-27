import xarray as xr
#import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
from datetime import timedelta
import gridpp
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
    

