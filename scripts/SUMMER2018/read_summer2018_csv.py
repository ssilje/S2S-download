
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




from S2S.file_handling import read_grib_file,read_grib_file_point, read_grib_file_slice_merge_ftype, check_file
from S2S.local_configuration import config

dates_monday = pd.date_range("20180503", periods=18, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20180507", periods=18, freq="7D") # forecasts start Thursday
dates_fcycle = dates_monday.union(dates_thursday)   


var_name_abbr='t2m' #tp

dirbase=config['csv_DIR']


file_hc =  '_'.join([dirbase]) + '/t2m_CY43R3_CY45R1_hindcast.csv'
file_fc =  '_'.join([dirbase]) + '/t2m_CY43R3_CY45R1_hindcast.csv'

print(file_hc)

hindcast_df = pd.read_csv(file_hc)

forecast_df = pd.read_csv(file_hc)

#for dd in dates_fcycle:
   # print(dd.strftime('%Y-%m-%d'))
#  hindcast[(hindcast["time"] ==d'd) & (hindcast["step"] =='7 days')]

d = dates_fcycle[0].strftime('%Y-%m-%d')

refyear = int(d[:4])
for i in range(refyear-20,refyear):
  hdate = d.replace('%i'%refyear,'%i'%i) 
  print(hdate)
  h_dd = hindcast[(hindcast["time"] ==hdate)]
 
 # h_dd[(h_dd["step"] =='7 days')] samme som h_dd[h_dd["step"].isin(['7 days','8 days'])]
  week2 = h_dd[h_dd["step"].isin(['7 days','8 days','9 days','10 days','11 days','12 days','13 days'])]

