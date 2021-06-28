import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import os
import scipy.stats as stats

from S2S.local_configuration import config
from S2S.data_handler        import BarentsWatch, ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.handle_datetime   as dt

class Hindcast:

    def __init__(
                    self,
                    var,
                    t_start,
                    t_end,
                    domainID,
                    high_res=False,
                    steps=None,
                    download=False,
                    process=False,
                    split_work=False,
                ):

        self.var            = var
        self.t_start        = t_start
        self.t_end          = t_end
        self.domainID       = domainID
        self.high_res       = high_res
        self.steps          = steps
        self.download       = download
        self.process        = process
        self.path           = config['VALID_DB']

        filename_absolute = self.filename_func('absolute')

        if self.process or not os.path.exists(self.path+filename_absolute):

            print('Process hindcast')
            if not os.path.exists(self.path):
                os.makedirs(self.path)

            if split_work:

                data_list = []

                t_end   = self.t_end
                t_start = self.t_start

                self.t_end = self.add_month(t_start)

                while self.smaller_than(self.t_end,t_end):

                    print('\tLoad hindcast')
                    raw = self.load_data()

                    print('\tApply 7D running mean along lead time dimension')
                    data = raw.rolling(step=7,center=True).mean()

                    if self.steps is not None:
                        print('\tKeep only specified lead times')
                        data = data.where(
                                        data.step.isin(self.steps),
                                        drop=True
                                    )

                    self.t_start = self.t_end
                    self.t_end   = self.add_month(self.t_end)
                    print(data)
                    data_list.append(data)

                self.data = xr.concat(data_list,'time')
                print(self.data)
            else:

                print('\tLoad hindcast')
                self.raw = self.load_data()

                print('\tApply 7D running mean along lead time dimension')
                self.data = self.raw.rolling(step=7,center=True).mean()

                if self.steps is not None:
                    print('\tKeep only specified lead times')
                    self.data = self.data.where(
                                            self.data.step.isin(self.steps),
                                            drop=True
                                        )

            self.data = self.drop_unwanted_dimensions(self.data)

            self.store(self.data,filename_absolute)

        self.data = self.load(filename_absolute)

        filename_anomalies = self.filename_func('anomalies')
        filename_mean      = self.filename_func('model_mean')
        filename_std       = self.filename_func('model_std')

        if self.process or not os.path.exists(self.path+filename_anomalies):

            print('\tCompute model climatology')
            self.mean,self.std = xh.c_climatology(self.data)

            self.mean = self.mean.rename(self.var)
            self.std  = self.std.rename(self.var)

            print('\tCompute anomalies')
            self.data_a = ( self.data - self.mean ) / self.std

            self.store(self.data_a,filename_anomalies)
            self.store(self.mean,filename_mean)
            self.store(self.std, filename_std)

        self.data_a = self.load(filename_anomalies)
        self.mean   = self.load(filename_mean)
        self.std    = self.load(filename_std)

    def load_data(self):

        return ECMWF_S2SH(high_res=self.high_res)\
                        .load(
                                self.var,
                                self.t_start,
                                self.t_end,
                                self.domainID,
                                self.download
                            )[self.var]-272.15

    @staticmethod
    def drop_unwanted_dimensions(data):
        try:
            data = data.drop('number')
        except ValueError:
            pass
        try:
            data = data.drop('surface')
        except ValueError:
            pass
        try:
            data = data.drop('valid_time')
        except ValueError:
            pass
        return data

    @staticmethod
    def add_month(time):

        current_year  = time[0]
        current_month = time[1]
        current_day   = time[2]

        if current_month%12==0:
            current_month = 1
            current_year += 1
        else:
            current_month += 1

        return (current_year,current_month,current_day)

    @staticmethod
    def smaller_than(low,high):
        if low[0]<high[0]:
            return True
        elif low[0]==high[0] and low[1]<=high[1]:
            return True
        else:
            return False

    def filename_func(self,filename):
        return '_'.join(
                        [
                            filename,
                            self.var,
                            dt.to_datetime(self.t_start).strftime('%Y-%m-%d'),
                            dt.to_datetime(self.t_end).strftime('%Y-%m-%d'),
                            self.domainID,
                            'hindcast'
                        ]
                    ) + '.nc'

    def store(self,file,filename):
        file.to_netcdf(self.path+filename)

    def load(self,filename):
        return xr.open_dataset(self.path+filename)[self.var]

class Observations:

    def __init__(
                    self,
                    name,
                    observations,
                    forecast,
                    process=False
                ):

        self.name           = name
        self.observations   = observations
        self.forecast       = forecast
        self.process        = process
        self.var            = forecast.data.name
        self.path           = config['VALID_DB']

        self.t_start        = (
                                observations.time.min().dt.year.values,
                                observations.time.min().dt.month.values,
                                observations.time.min().dt.day.values
                            )
        self.t_end          = (
                                observations.time.max().dt.year.values,
                                observations.time.max().dt.month.values,
                                observations.time.max().dt.day.values
                            )

        filename_absolute = self.filename_func('absolute')

        if self.process or not os.path.exists(self.path+filename_absolute):

            print('Process observations')
            if not os.path.exists(self.path):
                os.makedirs(self.path)

            print('\tAssign step dimension to observations')
            self.data = xh.at_validation(
                                        self.observations,
                                        forecast.data.time + forecast.data.step,
                                        ddays=1
                                        )

            self.data = self.data.rename(self.var)
            self.data = self.data.drop('validation_time')

            self.store(self.data,filename_absolute)

        self.data = self.load(filename_absolute)

        filename_anomalies = self.filename_func('anomalies')
        filename_mean      = self.filename_func('obs_mean')
        filename_std       = self.filename_func('obs_std')

        if self.process or not os.path.exists(self.path+filename_anomalies):

            print('\tCompute climatology')
            self.mean,self.std = xh.o_climatology(self.data)

            self.mean = self.mean.rename(self.var)
            self.std  = self.std.rename(self.var)

            self.data_a = ( self.data - self.mean ) / self.std

            self.store(self.data_a,filename_anomalies)
            self.store(self.mean,  filename_mean)
            self.store(self.std,   filename_std)

        self.data_a = self.load(filename_anomalies)
        self.mean   = self.load(filename_mean)
        self.std    = self.load(filename_std)

    def filename_func(self,filename):
        return '_'.join(
                        [
                            self.name,
                            filename,
                            dt.to_datetime(self.t_start).strftime('%Y-%m-%d'),
                            dt.to_datetime(self.t_end).strftime('%Y-%m-%d'),
                            'observations'
                        ]
                    ) + '.nc'

    def store(self,file,filename):
        file.to_netcdf(self.path+filename)

    def load(self,filename):
        return xr.open_dataset(self.path+filename)[self.var]

class Grid2Point:
    """
    observations and forecast must have common dimensions: time and step
    observations must be equipped with a location dimension
        (with arbitrary name)
    forecast must have member, lon and lat dimensions
    """

    def __init__(self,observations,forecast):

        self.observations = observations
        self.forecast     = forecast

    def correlation(self,step_dependent=False):
        return self.select_location_by_correlation(
                                                self.observations,
                                                self.forecast,
                                                step_dependent
                                            )

    def select_location_by_correlation(
                                        self,
                                        observations,
                                        forecast,
                                        step_dependent=False
                                    ):

        if step_dependent:
            input_core_dims  = [
                                ['time'],
                                ['member','time','lon','lat'],
                                ['member','time','lon','lat'],
                                ['time','lon','lat'],
                                ['time','lon','lat']
                               ]

            output_core_dims = [
                                    ['member','time'],
                                    ['member','time'],
                                    ['time'],
                                    ['time']
                                ]
        else:
            input_core_dims  = [
                                ['time','step'],
                                ['member','time','step','lon','lat'],
                                ['member','time','step','lon','lat'],
                                ['time','step','lon','lat'],
                                ['time','step','lon','lat']
                               ]

            output_core_dims = [
                                    ['member','time','step'],
                                    ['member','time','step'],
                                    ['time','step'],
                                    ['time','step']
                                ]

        o = observations.data_a

        try:
            o = o.drop(('lon','lat'))
        except ValueError:
            pass

        a,v,m,s,o = xr.align(
                                forecast.data_a,
                                forecast.data,
                                forecast.mean,
                                forecast.std,
                                o
                            )

        ds = xr.merge(
                        [
                            a.rename('a'),
                            v.rename('v'),
                            m.rename('m'),
                            s.rename('s'),
                            o.rename('o')
                    ],join='inner',compat='override'
                )

        a,v,m,s, = xr.apply_ufunc(
            self.get_highest_r, ds.o, ds.a, ds.v, ds.m, ds.s,
            input_core_dims  = input_core_dims,
            output_core_dims = output_core_dims,
            vectorize=True
            )

        forecast.data_a = a
        forecast.data   = v
        forecast.mean   = m
        forecast.std    = s

        return forecast

    @staticmethod
    def get_highest_r(o,a,v,m,s):
        """
        """
        # Flatten grid
        if len(a.shape) == 4: # step dependent
            a = a.reshape(a.shape[0],a.shape[1],-1)
            v = v.reshape(v.shape[0],v.shape[1],-1)
            m = m.reshape(m.shape[0],-1)
            s = s.reshape(s.shape[0],-1)

        else: # not step dependent
            a = a.reshape(a.shape[0],a.shape[1],a.shape[2],-1)
            v = v.reshape(v.shape[0],v.shape[1],v.shape[2],-1)
            m = m.reshape(m.shape[0],m.shape[1],-1)
            s = s.reshape(s.shape[0],s.shape[1],-1)

        # keep a with member dim
        A = a

        # mean over member dim
        a = a.mean(0)

        r = [] # for each gridpoint
        for ii in range(a.shape[-1]):

            # keep only nonnan forecast-observation pairs
            idx_bool = ~np.logical_or(
                                np.isnan(o),
                                np.isnan(a[...,ii])
                            )

            # if all nan (most likely landmask) give -99 as correlation coef.
            if idx_bool.sum()==0:
                r.append(-99)
            # otherwise pearsons r
            else:
                r.append(stats.pearsonr(o[idx_bool],a[idx_bool,ii])[0])

        highest_r = np.argmax(np.array(r))

        return  A[...,highest_r],v[...,highest_r],\
                m[...,highest_r],s[...,highest_r]
