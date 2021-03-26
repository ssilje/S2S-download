
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from calendar import monthrange,  monthcalendar, datetime
import gridpp
import json 
import os




from S2S.file_handling import read_grib_file,read_grib_file_point, read_grib_file_slice_merge_ftype, check_file
from S2S.local_configuration import config


var_name_abbr='t2m' #tp
#mdl_vrsn='CY43R3_CY45R1'
#S2S_dirbase=config['S2S_DIR_summer2018']
dirbase=config['csv_DIR']


file_hc =  '_'.join([dirbase]) + 't2m_CY43R3_CY45R1_hindcast.csv'

print(file_hc)

hindcast = pd.read_csv(file_hc)
