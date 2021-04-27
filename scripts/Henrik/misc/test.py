from .nird2local import pull
import xarray as xr
import numpy as np
import pandas as pd



# # TODO: make date range

# def date_range(start,end):
#     """
#     Returns list of dates between (and including) start and end.
#
#     args:
#         start: tuple of length 3
#         end:   tuple of length 3
#     """
#     for inst in :
#         for mm in range()

# sst_CY46R1_2020-05-18_hindcast.grb

def date_fmt(year,month,day):
    return '-'.join([str(year),str(month),str(day)])

var_name_abbr = 'sst'
model_version = 'CY46R1'
cast_type     = 'hindcast'

date  = date_fmt(
                    year='2020',
                    month='05',
                    day='18'
                )

filename = ['/' + '_'.join([var_name_abbr, model_version, date, cast_type]) + '.grb']
uni_path = '/SFE/hindcast/ECMWF/sst'
username = 'heau'

local_path = pull(uni_path,filename,username)
print(local_path+filename[0])

data = xr.open_dataset(
                        local_path+filename[0],
                        engine='cfgrib',
                        backend_kwargs={'filter_by_keys': {'dataType': 'cf'}}
                        )
print(data)

# print(data)
