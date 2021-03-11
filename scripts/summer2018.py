#%%
import sys  
sys.path.append('/nird/projects/NS9001K/sso102/Python/S2S-download/lib')  

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
import gridpp
import json 
import os

from globals import read_grib_file, read_grib_file_point, read_grib_file_merge_ftype

from settings_directories import DIR

#%% Dates
# var_name='sav300' 

lead_time=np.arange(1,47)
fcyear=2018
fcmonth=4
fcday=26

dates_monday = pd.date_range("20180426", periods=20, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20180430", periods=20, freq="7D") # forecasts start Thursday
dates_fcycle = dates_monday.union(dates_thursday)   


#dates_fcycle=pd.date_range(start=f'{fcyear}-{fcmonth}-{fcday}', periods=2, freq='7D') # forecasts start Monday

#%% Read in data for a given date

var_name_abbr='tp'
mdl_vrsn='CY43R3_CY45R1''
S2S_dirbase=DIR['S2SEUR_DIR']
#product='forecast'
curr_date=dates_fcycle[1].strftime('%Y-%m-%d')

for product in (
        'forecast',
        'hindcast',
    ):
        name = '%s_%s'%('ds',product)
        name = read_grib_file_merge_ftype(
        S2S_dirbase=S2S_dirbase,
        product=product,
        model_version=mdl_vrsn,
        var_name_abbr=var_name_abbr,
        date_str=curr_date,
        )
        print(name.head())

