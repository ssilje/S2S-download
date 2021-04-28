# environment dependencies
import numpy  as np
import xarray as xr
import pandas as pd

# local dependencies
from S2S.local_configuration import config
from .handle_domain import get_bounds
import scripts.Henrik.handle_datetime as dt

class SST:
    """
    Loads daily: S2S data, hindcast and forecast, and ERA5.
    """

    def __init__(self,prnt=True):

        self.ERA5   = xr.Dataset()
        self.S2SH   = xr.Dataset()
        self.S2SF   = xr.Dataset()

        self.labels = [
                        'ERA5',
                        'S2SH',
                        'S2SF'
                        ]

        self.ftype = {
                        'ERA5':'reanalysis',
                        'S2SH':'hindcast',
                        'S2SF':'forecast'
                        }

        self.loaded = dict.fromkeys(self.labels)
        self.path   = dict.fromkeys(self.labels)
        self.time   = dict.fromkeys(self.labels)

        self.domainID = 'full_grid'
        self.bounds   = ()

    def load(
                self,
                data_label,
                t_start,
                t_end,
                domainID=None,
                match_fc=False,
                fc_per_week=1
            ):

        self.domainID = domainID if domainID else self.domainID
        self.bounds   = get_bounds(domainID)

        # if not self.is_loaded(data_label) or reload:

        self.execute_loading_sequence(
                                        data_label,
                                        t_start,
                                        t_end,
                                        match_fc,
                                        fc_per_week
                                    )

    def execute_loading_sequence(
                                    self,
                                    data_label,
                                    t_start,
                                    t_end,
                                    match_fc=False,
                                    fc_per_week=1
                                ):

        self.path[data_label] = config[data_label]
        bounds                = self.bounds

        if data_label=='ERA5':

            self.time[data_label] = dt.days_from(t_start,t_end,match_fc)

            chunk = []

            for date in self.time[data_label]:

                chunk.append(
                    self.load_era(
                        self.path[data_label],
                        date,
                        self.bounds
                        )
                    )

            self.ERA5 = xr.concat(chunk,'time') # existing dimension

        else:

            if fc_per_week==1:
                self.time[data_label] = dt.weekly_forecast_cycle(t_start,t_end)
            else:
                self.time[data_label] = dt.forecast_cycle(t_start,t_end)

            chunk = []

            for date in self.time[data_label]:

                chunk.append(
                        self.load_S2S(
                            self.path[data_label],
                            date,
                            self.bounds,
                            self.ftype[data_label]
                            )
                        )

            if data_label=='S2SH':

                self.S2SH = xr.concat(chunk,'time') # existing dimension

            else:

                self.S2SF = xr.concat(chunk,'cast_time') # new dimension

    def is_loaded(self,data_label):
        if self.loaded[data_label]:
            return True
        else:
            return False

    @staticmethod
    def load_era(path,date,bounds):

        filename  = '_'.join(['sea_surface_temperature',date.strftime('%Y%m%d')])+'.nc'

        open_data = xr.open_dataset(path+filename)
        open_data = open_data.sel(
                        lat=slice(bounds[2],bounds[3]),
                        lon=slice(bounds[0],bounds[1])
                        )
        open_data = open_data.resample(time='D').mean()

        return open_data

    @staticmethod
    def load_S2S(path,date,bounds,ftype):

        filename = '_'.join(['sst','CY46R1_CY47R1',date.strftime('%Y-%m-%d'),'pf',ftype])+'.grb'
        open_data = xr.open_dataset(path+filename,engine='cfgrib')
        open_data = open_data.sortby('latitude', ascending=True)
        open_data = open_data.sel(
                        latitude=slice(bounds[2],bounds[3]),
                        longitude=slice(bounds[0],bounds[1])
                        )

        return open_data
