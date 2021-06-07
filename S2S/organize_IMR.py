import pandas as pd
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt

from S2S.local_configuration import config

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
stations  = [
                'Bud',
                'Eggum',
                'Ingoy',
                'Indre Utsira',
                'Lista',
                'Skrova',
                'Sognesjoen',
                'Ytre Utsira'
            ]

coords = dict.fromkeys(stations)
for n,filename in enumerate(filenames):

    in_data = xr.open_dataset(config['IMR']+filename+'.nc')
    lon     = in_data.LONGITUDE.values[0]
    lat     = in_data.LATITUDE.values[0]

    coords[stations[n]] = {'label':filename,'lon':lon,'lat':lat}

data  = dict.fromkeys(stations)
for key in data:
    data[key] = {
                'sst':[],
                'time':[],
                'lon':[],
                'lat':[]
                }
with open(config['IMR']+'stasjonsdata.txt','rb') as file:
    for n,line in enumerate(file):
        if n>0:

            try:
                line = line.decode()
                line = line.split()

                if len(line)==9:
                    station = line.pop(0)+' '+line.pop(0)
                else:
                    station = line.pop(0)

                try:
                    coords[station]

                    date    = line[0]
                    time    = line[1]
                    temp    = float(line[3])

                    time = pd.Timestamp(date + ' ' + time)

                    data[station]['sst'].append(temp)
                    data[station]['time'].append(time)

                except KeyError:
                    pass

            except UnicodeDecodeError:
                pass

for station in data:
    ds = xr.DataArray(
                        data=data[station]['sst'],
                        dims=['time'],
                        coords=dict(
                                location=coords[station]['label'],
                                time=data[station]['time'],
                                lon=coords[station]['lon'],
                                lat=coords[station]['lat']
                                )
                        ).rename('sst')
    ds.to_netcdf(config['IMR']+coords[station]['label']+'_organized.nc')
