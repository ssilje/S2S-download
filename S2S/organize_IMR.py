import pandas as pd
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt

from S2S.local_configuration import config

local_path = config['IMR']
filenames  = [
                'bud',
                'eggum',
                'ingoy',
                'iutsira',
                'lista',
                'skrova',
                'sognesjoen',
                'yutsira'
            ]

start = pd.Timestamp('01-01-2000')
print('With values from',start,'\n')
out = []
for filename in filenames:
    data = xr.open_dataset(local_path+filename+'.nc')

    for d in [0]:

        data = data.sel(
                DEPTH=d,
                TIME=slice(start,data.TIME.isel(TIME=-1))
                )

        lon = data.LONGITUDE.values[0]
        lat = data.LATITUDE.values[0]

        # Change to adjusted?
        data = data.TEMP.where(data.TEMP_QC<3)
        # time = data.TIME.where(data.TIME_QC<3)

        data = data.rename({'TIME':'time'}).rename('sst')
        data = data.assign_coords(location=filename)
        data = data.assign_coords(lon=lon)
        data = data.assign_coords(lat=lat)

    out.append(data)

    print('\tProcessing',filename,'-\t',len(data.values),'observations')

xr.concat(out,'location').to_netcdf(local_path+'fastestasjoner.nc')
print('\n...to',local_path+'fastestasjoner.nc')
