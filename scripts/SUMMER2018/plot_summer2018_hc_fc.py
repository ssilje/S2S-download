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
    


ds_hindcast_stat= ds_hindcast.groupby(["step"])[["t2m"]].describe()
ds_forecast_stat= ds_forecast.groupby(["step"])[["t2m"]].describe()

for w in fc_week:
    temp_hc = ds_hindcast_stat.loc[fc_week[w]].mean()
    temp_hc_df=pd.DataFrame([temp_hc],index=[w])
     
    temp_fc = ds_forecast_stat.loc[fc_week[w]].mean()
    temp_fc_df=pd.DataFrame([temp_fc],index=[w])

    if w == "week1":
        hindcast_stat = temp_hc_df
        forecast_stat = temp_fc_df
      
    else:
        hindcast_stat = hindcast_stat.append(temp_hc_df)
        forecast_stat = forecast_stat.append(temp_fc_df)
       
 
stats_hc = [{
    "label": 't2m-hc week1',  # not required
    "med": hindcast_stat.t2m["mean"]["week1"], #5.5
    "q1": hindcast_stat.t2m["25%"]["week1"],
    "q3": hindcast_stat.t2m["75%"]["week1"],
    # "cilo": 5.3 # not required
    # "cihi": 5.7 # not required
    "whislo": hindcast_stat.t2m["min"]["week1"],  # required
    "whishi": hindcast_stat.t2m["max"]["week1"],  # required
    "fliers": []  # required if showfliers=True
    }]

stats_fc = [{
    "label": 't2m-fc week1',  # not required
    "med": forecast_stat.t2m["mean"]["week1"], #5.5
    "q1": forecast_stat.t2m["25%"]["week1"],
    "q3": forecast_stat.t2m["75%"]["week1"],
    # "cilo": 5.3 # not required
    # "cihi": 5.7 # not required
    "whislo": forecast_stat.t2m["min"]["week1"],  # required
    "whishi": forecast_stat.t2m["max"]["week1"],  # required
    "fliers": []  # required if showfliers=True
    }]

fs = 10  # fontsize

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
axes.bxp(1,stats_hc,patch_artist=True,)
axes.bxp(2,stats_fc)
axes.set_title('Boxplot for precalculated statistics', fontsize=fs)
fig.savefig('test.png')

#test.mean
##test.t2m('mean')
#test.t2m['mean']
#test.t2m['25%']


#d = curr_date

#refyear = int(d[:4])
#for i in range(refyear-20,refyear):
 # hdate = d.replace('%i'%refyear,'%i'%i) 
 # print(hdate)
 # h_dd = hindcast[(hindcast["time"] ==hdate)]
 
 # h_dd[(h_dd["step"] =='7 days')] samme som h_dd[h_dd["step"].isin(['7 days','8 days'])]
 # week2 = h_dd[h_dd["step"].isin(['7 days','8 days','9 days','10 days','11 days','12 days','13 days'])]
