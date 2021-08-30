import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import os
import scipy.stats as stats

from S2S.local_configuration import config
from S2S.data_handler        import BarentsWatch, ERA5, ECMWF_S2SH, ECMWF_S2SF, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.handle_datetime   as dt

class Hindcast:
    """
    Loads hindcast from S2S database, computes weekly means and then provides

        self.data:      the absolute values of the hindcast (xarray.DataArray)
        self.data_a:    the anomlies of the data relative to model climatology
                        (xarray.DataArray)
        self.mean:      the model mean climatology (over 30-day running window)
                        (xarray.DataArray)
        self.std:       the model std climatology (over 30-day running window)
                        (xarray.DataArray)

    Arguments to __init__

        var:        the name of the variable, should correspond to filenames in
                    the S2S database (string)
        t_start:    start time of files to load (tuple of int; (year,monty,day))
        t_end:      end time of files to load (tuple of int; (year,monty,day))
                    both start and end are included to load files
        bounds:     bounds of lon-lat grid to load
                    (tuple of float; (min lon, max lon, min lat, max lat))
                    Current grid is 0-360 in longitude direction starting at 0
                    at the exact location of Boris Johnson (does not vary much).
        high_res:   changes path to 0.5 degree grid hindcast (only for SST).
                    Default is 1.5 degree grid.
        steps:      The steps one would like to store (after computing running
                    7-day means along step (lead time) dimension).
                    In combination with split_work=True, this could save some
                    internal memory, by re-arrangeing(?) the work flow and
                    tossing out lead times at an earlier stage. Default keeps
                    all lead times.
        dowload:    if True, and process is set to True, forces download from
                    S2S database even if temporary files are found. Default is
                    False.
        process:    if True, forces the computation of model climatology and
                    anomalies. Default is False.
        split_work: if True, exploits the steps given in argument 'steps' to
                    reduce the use of internal memory. Default is False.
        period:     The period that climatology should be computed over
                    (list of tuples of int; [(year,monty,day),(year,monty,day)])
                    Useful if the the observational period is shorter than the
                    hindcast period. Default is None.
    """
    def __init__(
                    self,
                    var,
                    t_start,
                    t_end,
                    bounds,
                    high_res   = False,
                    steps      = None,
                    download   = False,
                    process    = False,
                    split_work = False,
                    period     = None,
                    cross_val  = False
                ):

        self.var            = var
        self.t_start        = t_start
        self.t_end          = t_end
        self.bounds         = bounds
        self.high_res       = high_res
        self.steps          = steps
        self.download       = download
        self.process        = process
        self.path           = config['VALID_DB']
        self.period         = period
        self.cross_val      = cross_val

        filename_absolute = self.filename_func('absolute')

        if self.process or not os.path.exists(self.path+filename_absolute):

            print('Process hindcast')
            if not os.path.exists(self.path):
                os.makedirs(self.path)

            if split_work:

                data_list = []

                # store the actual start and end times
                t_end   = self.t_end
                t_start = self.t_start

                # assign new end time, one month after start time
                self.t_end = self.add_month(t_start)

                # until end time is reached; load one month at the time
                while self.smaller_than(self.t_end,t_end):

                    print('\tLoad hindcast')
                    raw = self.load_data().sortby('time')

                    # Select only given period
                    if self.period is not None:

                        self.raw = self.raw.sel(
                                        time=slice(
                                                dt.to_datetime(self.period[0]),
                                                dt.to_datetime(self.period[1])
                                                )
                                            )

                    print('\tApply 7D running mean along lead time dimension')
                    data = raw.rolling(step=7,center=True).mean()

                    if self.steps is not None:
                        print('\tKeep only specified lead times')
                        data = data.where(
                                        data.step.isin(self.steps),
                                        drop=True
                                    )

                    # update times to the next month
                    self.t_start = self.t_end
                    self.t_end   = self.add_month(self.t_end)

                    # append to list
                    data_list.append(data)

                # concatinate xarray.DataArrays in list to one xarray.DataArray
                self.data = xr.concat(data_list,'time')

                # deal with duplicates along time dimesion
                self.data = self.data.groupby('time').mean()

                # restore original times of loading
                self.t_start = t_start
                self.t_end   = t_end

            else:

                print('\tLoad hindcast')
                self.raw = self.load_data().sortby('time')

                # Select only given period
                if self.period is not None:

                    self.raw = self.raw.sel(
                                    time=slice(
                                            dt.to_datetime(self.period[0]),
                                            dt.to_datetime(self.period[1])
                                            )
                                        )

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
            self.mean,self.std = xh.c_climatology(
                                            self.data,
                                            cross_validation = self.cross_val
                                                )

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

        data = ECMWF_S2SH(high_res=self.high_res)\
                        .load(
                                self.var,
                                self.t_start,
                                self.t_end,
                                self.bounds,
                                self.download
                            )[self.var]

        # Converts Kelvin to degC, if this should be done here can be discussed?
        if self.var == 'sst':
            data = data - 273.15

        return data

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
                            '%s_%s-%s_%s'%(self.bounds),
                            'hindcast'
                        ]
                    ) + '.nc'

    def store(self,file,filename):
        file.to_netcdf(self.path+filename)

    def load(self,filename):
        return xr.open_dataset(self.path+filename)[self.var]

class Forecast:
    """
    Loads forecasts from S2S database, computes weekly means and then provides

        self.data:      the absolute values of the forecast (xarray.DataArray)
        self.data_a:    the anomlies of the data relative to model hindcast climatology
                        (xarray.DataArray)
        self.hc_mean:   the model hindcast mean climatology (over 30-day running window)
                        (xarray.DataArray)
        self.hc_std:    the model hindcast std climatology (over 30-day running window)
                        (xarray.DataArray)
        self.hc_abs:    absolute values of the 20 years of hindcasts corresponding to the
                        forecast for the given initialization(s) (xarray.DataArray)
        self.hc_anom:   hindcast anomalies relative to the same model hindcast climatology
                        as used for forecast anomalies (xarray.DataArray)

    Arguments to __init__

        var:        the name of the variable, should correspond to filenames in
                    the S2S database (string)
        t_start:    start time of files to load (tuple of int; (year,monty,day))
        t_end:      end time of files to load (tuple of int; (year,monty,day))
                    both start and end are included to load files
        bounds:     bounds of lon-lat grid to load
                    (tuple of float; (min lon, max lon, min lat, max lat))
                    Current grid is 0-360 in longitude direction starting at 0
                    at the exact location of Boris Johnson (does not vary much).
        high_res:   changes path to 0.5 degree grid hindcast (only for SST).
                    Default is 1.5 degree grid.
                    (high_res not currently available for forecasts ---2021-07-26---)
        steps:      The steps one would like to store (after computing running
                    7-day means along step (lead time) dimension).
                    In combination with split_work=True, this could save some
                    internal memory, by re-arrangeing(?) the work flow and
                    tossing out lead times at an earlier stage. Default keeps
                    all lead times.
        dowload:    if True, and process is set to True, forces download from
                    S2S database even if temporary files are found. Default is
                    False.
        process:    if True, forces the computation of model climatology and
                    anomalies. Default is False.
        split_work: if True, exploits the steps given in argument 'steps' to
                    reduce the use of internal memory. Default is False.
    """
    def __init__(
                    self,
                    var,
                    t_start,
                    t_end,
                    bounds,
                    high_res=False,
                    steps=None,
                    download=False,
                    process=False,
                    split_work=False,
                ):

        self.var            = var
        self.t_start        = t_start
        self.t_end          = t_end
        self.bounds         = bounds
        self.high_res       = high_res
        self.steps          = steps
        self.download       = download
        self.process        = process
        self.path           = config['VALID_DB']

        filename_absolute = self.filename_func('absolute')

        if self.process or not os.path.exists(self.path+filename_absolute):

            print('Process forecast')
            if not os.path.exists(self.path):
                os.makedirs(self.path)

            if split_work:

                data_list = []

                # store the actual start and end times
                t_end   = self.t_end
                t_start = self.t_start

                # assign new end time, one month after start time
                self.t_end = self.add_month(t_start)

                # until end time is reached; load one month at the time
                while self.smaller_than(self.t_end,t_end):

                    print('\tLoad forecast')
                    raw = self.load_data()

                    print('\tApply 7D running mean along lead time dimension')
                    data = raw.rolling(step=7,center=True).mean()

                    if self.steps is not None:
                        print('\tKeep only specified lead times')
                        data = data.where(
                                        data.step.isin(self.steps),
                                        drop=True
                                    )

                    # update times to the next month
                    self.t_start = self.t_end
                    self.t_end   = self.add_month(self.t_end)

                    # append to list
                    data_list.append(data)

                # concatinate xarray.DataArrays in list to one xarray.DataArray
                self.data = xr.concat(data_list,'time')

                # # deal with duplicates along time dimesion
                # self.data = self.data.groupby('time').mean()

                # restore original times of loading
                self.t_start = t_start
                self.t_end   = t_end

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
        filename_anomalies_nonstd = self.filename_func('anomalies_nonstd')
        filename_mean      = self.filename_func('model_mean')
        filename_std       = self.filename_func('model_std')

        if self.process or not os.path.exists(self.path+filename_anomalies):

            # load the hindcasts in order to get the model climatology:
            print('\tloading corresponding hindcasts for climatology')
            hc_ref = Hindcast(
                                var                 = self.var,
                                t_start             = self.t_start,
                                t_end               = self.t_end,
                                bounds              = self.bounds,
                                high_res            = self.high_res,
                                steps               = self.steps,
                                download            = self.download,
                                process             = self.process,
                                split_work          = split_work
            )

            # climatologies saved in hc_ref.
            # subtract clim mean and divide by clim std.
            # climatology files are extended to a value each year are (should be) the same for every year
            # so it should be possible to take the mean over all hindcast years (try to find better solution):
            self.data_anom = (self.data.groupby('time.dayofyear') - hc_ref.mean.groupby('time.dayofyear').mean())
            self.data_a = self.data_anom.groupby('time.dayofyear')/hc_ref.std.groupby('time.dayofyear').mean()

            self.store(self.data_a,filename_anomalies)
            self.store(self.data_anom,filename_anomalies_nonstd)
            self.store(hc_ref.mean,filename_mean)
            self.store(hc_ref.std,filename_std)

        self.data_a = self.load(filename_anomalies)
        self.data_anom   = self.load(filename_anomalies_nonstd)
        self.hc_mean = self.load(filename_mean)
        self.hc_std = self.load(filename_std)
        self.hc_abs = hc_ref.data
        self.hc_anom = hc_ref.data_a

    def load_data(self):

        data = ECMWF_S2SF(high_res=self.high_res)\
                        .load(
                                self.var,
                                self.t_start,
                                self.t_end,
                                self.bounds,
                                self.download
                            )[self.var]

        # Converts Kelvin to degC, if this should be done here can be discussed?
        if self.var == 'sst':
            data = data - 273.15

        return data

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
                            '%s_%s-%s_%s'%(self.bounds),
                            'forecast'
                        ]
                    ) + '.nc'

    def store(self,file,filename):
        file.to_netcdf(self.path+filename)

    def load(self,filename):
        return xr.open_dataset(self.path+filename)[self.var]

class Observations:
    """
    Stacks observations along step (lead time) dimension to match forecast
    for matrix operations. Further computes climatology using a 30-day running
    window with cross validation (leaving out the current year). Produces:

        self.data:      the absolute values of the observations
                        (xarray.DataArray)
        self.data_a:    the anomalies of the data relative to climatology
                        (xarray.DataArray)
        self.mean:      the mean climatology (over 30-day running window)
                        (xarray.DataArray)
        self.std:       the std climatology (over 30-day running window)
                        (xarray.DataArray)
        self.init_a     the observed value at forecast initialization time,
                        stacked along step (lead time) dimension. Given as
                        anomlies.

    Arguments to __init__

        name:           name of the observation data, to create temporary files
        observations:   xarray.DataArray with required dimension time. Must
                        be seven day means.
        forecast:       process.Hindcast like object. Assimilates the
                        observations to the attributes of this input.
        process:        Forces data processing even if temporary files are
                        found. Temporary file names are not so flexible at the
                        moment, can be smart to keep this option on if input
                        observations are changed.
    """
    def __init__(
                    self,
                    name,
                    observations,
                    forecast,
                    process=False
                ):

        self.name           = name
        self.observations   = observations.sortby('time')
        self.forecast       = forecast
        self.process        = process
        self.var            = forecast.data.name
        self.path           = config['VALID_DB']

        self.forecast.data  = self.forecast.data.sortby(['time','step'])

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

            # lon/lat are lost in storage process
            # if contained in encoding['coordinates']
            coords = self.data.encoding['coordinates'].split(' ')
            while 'lon' in coords: coords.remove('lon')
            while 'lat' in coords: coords.remove('lat')
            self.data.encoding['coordinates'] = ' '.join(coords)

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

        filename_init = self.filename_func('init')

        if self.process or not os.path.exists(self.path+filename_init):

            print('\tGather observations at model initalization')
            init_mean,init_std = xh.o_climatology(self.observations)

            init_mean = init_mean.rename(self.var)
            init_std  = init_std.rename(self.var)

            ####################################################################
            # Deals with duplicates along time dimensions (can occur in
            # stacking - restacking occuring o_climatology, does occur for ERA5)
            _,c = np.unique(init_mean.time.values, return_counts=True)
            if len(c[c>1])>0:
                init_mean = init_mean.groupby('time').mean(skipna=True)

            _,c = np.unique(init_std.time.values, return_counts=True)
            if len(c[c>1])>0:
                init_std = init_std.groupby('time').mean(skipna=True)
            ####################################################################

            init_obs_a = ( self.observations - init_mean ) / init_std

            init_obs_a = init_obs_a.reindex(
                            {'time':forecast.data.time},
                            method='pad',
                            tolerance='7D'
                        ).broadcast_like(self.data)

            self.store(init_obs_a,filename_init)

        self.init_a = self.load(filename_init)

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
    Class rountine correlation() returns a forecast (process.Hindcast like)
    with corresponding dimensions to observations.
    The routine selects the highest correlated (Pearson) gridpoint to represent
    each dimension in observations, respectively.

    Arguments to __init__

        observations:   process.Observations like object
        forecast:       process.Hindcast like object

        Note:
            - observations and forecast must have common dimensions:
            time and step
            - forecast must have member, lon and lat dimensions
    """

    def __init__(self,observations,forecast):

        self.observations = observations
        self.forecast     = forecast

    def correlation(self,step_dependent=False):
        """
        Returns a forecast (process.Hindcast like) with corresponding dimensions
        to observations. The routine selects the highest correlated (Pearson)
        gridpoint to represent each dimension in observations, respectively.

        args:
            step_dependent: if True, correlation is done respectively of each
                            step (lead time). Default is False.

        returns:
            process.Hindcast like object with 'regridded' to dimensions of
            observations, additional to time and step.
        """
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
            idx_bool = np.logical_and(
                                np.isfinite(o),
                                np.isfinite(a[...,ii])
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
