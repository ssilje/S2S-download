import os
import numpy as np
import xarray as xr

from S2S.local_configuration import config
import scripts.Henrik.handle_datetime as dt
from .create_domain_file import make_dir

from .data_handler import SST

def get_observations(domainID,start,end,download=0):

    obspath = '%s%s_%s-%s_%s_%s%s'%(
                            config['VALID_DB'],
                            'sst',
                            dt.to_datetime(start).strftime('%Y-%m-%d'),
                            dt.to_datetime(end).strftime('%Y-%m-%d'),
                            'reanalysis',
                            domainID,
                            '.nc'
                            )

    if not download and not os.path.exists(obspath):
        download = 1

    while True:

        if download:

            observations = SST()
            observations.load('ERA5',start,end,domainID)

            make_dir('/'.join(obspath.split('/')[:-1]))
            observations.ERA5.transpose('time','lon','lat').to_netcdf(obspath)

            download = 0

        else:

            return xr.open_dataset(obspath)

def get_hindcast(domainID,start,end,download=0):

    hc_filepath = '%s%s_%s-%s_%s_%s%s'%(
                            config['VALID_DB'],
                            'sst',
                            dt.to_datetime(start).strftime('%Y-%m-%d'),
                            dt.to_datetime(end).strftime('%Y-%m-%d'),
                            'hc',
                            domainID,
                            '.nc'
                            )

    if not download and not os.path.exists(hc_filepath):
        download = 1

    while True:

        if download:

            sst = SST()
            sst.load('S2SH',start,end,domainID)

            hindcast = sst.S2SH.transpose(
                                        'number','step','time',
                                        'longitude','latitude'
                                        )

            hindcast = hindcast.rename(
                                        {
                                            'number':'member',
                                            'longitude':'lon',
                                            'latitude':'lat'
                                            }
                                        )

            make_dir('/'.join(hc_filepath.split('/')[:-1]))
            hindcast.to_netcdf(hc_filepath)

            download = 0

        else:

            return xr.open_dataset(hc_filepath)
