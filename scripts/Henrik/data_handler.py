# environment dependencies
import numpy  as np
import xarray as xr
import pandas as pd
import os

# local dependencies
from S2S.local_configuration import config
from .handle_domain import get_bounds, get_domain
import scripts.Henrik.handle_datetime as dt
import scripts.Henrik.organize_barentzwatch as organize_barentzwatch

class Archive:
    """
    Returns local paths and filename functions
    """

    def __init__(self):

        self.source = ['ERA5','S2SH','S2SF','BW']
        self.path   = dict((source,config[source]) for source in self.source)
        self.ftype  = {
                        'ERA5':'reanalysis',
                        'S2SH':'hindcast',
                        'S2SF':'forecast'
                        }

        self.dimension_parser = {
                        'number':'member',
                        'longitude':'lon',
                        'latitude':'lat'
                        }

    @staticmethod
    def ERA5_in_filename(**kwargs):

        var  = kwargs['var']
        date = kwargs['date']

        var = {
                'sst':'sea_surface_temperature'
            }[var]
        return '_'.join([var,date.strftime('%Y%m%d')])+'.nc'

    @staticmethod
    def S2S_in_filename(**kwargs):

        var   = kwargs['var']
        date  = kwargs['date']
        run   = kwargs['run']
        ftype = kwargs['ftype']

        if kwargs['high_res']:
            return '_'.join(
                    [
                        var,
                        'CY46R1_CY47R1',
                        '05x05',
                        date.strftime('%Y-%m-%d'),
                        run,
                        ftype
                    ]
                ) + '.grb'
        else:
            return '_'.join(
                    [
                        var,
                        'CY46R1_CY47R1',
                        date.strftime('%Y-%m-%d'),
                        run,
                        ftype
                    ]
                ) + '.grb'

    @staticmethod
    def BW_in_filename(**kwargs):

        location = kwargs['location']

        return '_'.join(['barentswatch',str(location)])+'.nc'

    @staticmethod
    def out_filename(var,start,end,domainID,label):
        return '%s_%s-%s_%s_%s%s'%(
                            var,
                            dt.to_datetime(start).strftime('%Y-%m-%d'),
                            dt.to_datetime(end).strftime('%Y-%m-%d'),
                            label,
                            domainID,
                            '.nc'
                        )

    @staticmethod
    def make_dir(path):
        """
        Creates directory if it does not exist.

        args:
            path: str
        """
        if not os.path.exists(path):
            os.makedirs(path)

class LoadLocal:
    """
    Base class for loading local data
    """

    def __init__(self):

        self.prnt            = True

        self.label           = ''
        self.var             = ''
        self.download        = False
        self.loading_options = {
                                'label':            None,
                                'load_time':        None,
                                'concat_dimension': None,
                                'resample':         None,
                                'sort_by':          None,
                                'control_run':      None,
                                'engine':           None,
                                'high_res':         None
                            }

        self.out_path     = ''
        self.in_path      = ''

        self.out_filename = ''

        self.start_time   = None
        self.end_time     = None

        self.domainID     = 'full_grid'
        self.bounds       = ()

        self.dimension_parser = Archive().dimension_parser

    def rename_dimensions(self,ds):

        for dim in ds.dims:
            try:
                ds = ds.rename({dim:self.dimension_parser[dim]})
            except KeyError:
                pass
        return ds

    def filename(self,key='in'):

        label = self.label

        if label=='ERA5' and key=='in':
            return Archive().ERA5_in_filename

        if label=='S2SH' and key=='in':
            return Archive().S2S_in_filename

        if label=='S2SF' and key=='in':
            return Archive().S2S_in_filename

        if label=='BW' and key=='in':
            return Archive().BW_in_filename

        if key=='out':
            return Archive().out_filename

    def load_frequency(self):

        option = self.loading_options['load_time']

        if option=='daily':
            return dt.days_from(self.start_time,self.end_time)

        elif option=='weekly_forecast_cycle':
            return dt.weekly_forecast_cycle(self.start_time,self.end_time)

        elif option=='forecast_cycle':
            return dt.forecast_cycle(self.start_time,self.end_time)
        else:
            raise KeyError
            exit()

    def execute_loading_sequence(self):

        chunk = []

        sort_by     = self.loading_options['sort_by']
        resample    = self.loading_options['resample']
        engine      = self.loading_options['engine']
        dimension   = self.loading_options['concat_dimension']
        control_run = self.loading_options['control_run']
        high_res    = self.loading_options['high_res']
        ftype       = Archive().ftype[self.label]

        filename_func = self.filename(key='in')

        for time in self.load_frequency():

            runs   = ['pf','cf'] if control_run else [None]

            OK = True
            for n,run in enumerate(runs):

                filename = filename_func(
                                        var      = self.var,
                                        date     = time,
                                        run      = run,
                                        ftype    = ftype,
                                        high_res = high_res
                                    )

                OK = os.path.exists(self.in_path+filename)
            
            if OK:

                members = []
                for n,run in enumerate(runs):

                    filename = filename_func(
                                            var      = self.var,
                                            date     = time,
                                            run      = run,
                                            ftype    = ftype,
                                            high_res = high_res
                                        )

                    if n>0:
                        members.append(open_data)


                    open_data = xr.open_dataset(self.in_path+\
                                                    filename,engine=engine)

                    open_data = self.rename_dimensions(open_data)

                    if sort_by:
                        open_data = open_data.sortby(sort_by,ascending=True)

                    open_data = open_data.sel(
                                    lat=slice(self.bounds[2],self.bounds[3]),
                                    lon=slice(self.bounds[0],self.bounds[1])
                                    )

                    if resample:
                        open_data = open_data.resample(time=resample).mean()

                    if run=='cf':
                        open_data = open_data.expand_dims('member')\
                                        .assign_coords(member=pd.Index([0]))

                    if self.prnt:
                        print(filename)

                if n>0:
                    members.append(open_data)
                    open_data = xr.concat(members,'member')

                chunk.append(open_data)

        return xr.concat(chunk,dimension)

    def load(self,var,start_time,end_time,domainID,download=False,prnt=True):

        archive = Archive()

        self.prnt         = prnt

        self.var          = var
        self.start_time   = start_time
        self.end_time     = end_time
        self.domainID     = domainID
        self.bounds       = get_bounds(domainID)
        self.download     = download
        self.out_filename = archive.out_filename(
                                                var=var,
                                                start=start_time,
                                                end=end_time,
                                                domainID=domainID,
                                                label=archive.ftype[self.label]
                                                )

        while not self.download and not os.path.exists(self.out_path
                                                        +self.out_filename):

            archive.make_dir(self.out_path)

            data = self.execute_loading_sequence()

            if self.label=='ERA5':
                data.transpose('time','lon','lat').to_netcdf(self.out_path
                                                            +self.out_filename)
            elif self.label=='S2SH':
                data.transpose(
                            'member','step','time',
                            'lon','lat'
                            ).to_netcdf(self.out_path+self.out_filename)
            else: # TODO make option for forecast

                data.to_netcdf(self.out_path+self.out_filename)

            self.download = False

        return xr.open_dataset(self.out_path+self.out_filename)

class ERA5(LoadLocal):

    def __init__(self,high_res=False):

        super().__init__()

        self.label           = 'ERA5'

        self.loading_options['load_time']        = 'daily'
        self.loading_options['concat_dimension'] = 'time'
        self.loading_options['resample']         = 'D'
        self.loading_options['sort_by']          = 'lat'
        self.loading_options['control_run']      = False
        self.loading_options['engine']           = 'netcdf4'

        if high_res:
            self.loading_options['high_res']     = True
            self.in_path                         = config[self.label+'_HR']
        else:
            self.in_path                         = config[self.label]

        self.out_path                            = config['VALID_DB']

class ECMWF_S2SH(LoadLocal):

    def __init__(self,high_res=False):

        super().__init__()

        self.label           = 'S2SH'

        self.loading_options['load_time']        = 'weekly_forecast_cycle'
        self.loading_options['concat_dimension'] = 'time'
        self.loading_options['resample']         = False
        self.loading_options['sort_by']          = 'lat'
        self.loading_options['control_run']      = True
        self.loading_options['engine']           = 'cfgrib'

        if high_res:
            self.loading_options['high_res']     = True
            self.in_path                         = config[self.label+'_HR']
        else:
            self.in_path                         = config[self.label]

        self.out_path                            = config['VALID_DB']


class BarentsWatch:

    def __init__(self,prnt=True):

        self.labels = [
                        'DATA',
                        ]

        self.ftype = {
                        'DATA':'point_measurement',
                        }

        self.loaded = dict.fromkeys(self.labels)
        self.path   = dict.fromkeys(self.labels)
        self.time   = dict.fromkeys(self.labels)

        self.path['DATA'] = config['BW']

    def load(self,location,data_label='DATA'):

        archive = Archive()

        if not hasattr(location,'__iter__'):
            if isinstance(location,str) or isinstance(location,int):
                location = [location]
            else:
                raise ValueError

        location_2 = []
        if isinstance(location[0],str):
            for loc in location:
                location_2.append(get_domain(loc)['localityNo'])

        location = location_2

        if not os.path.exists(self.path['DATA']+\
                                archive.BW_in_filename(location=location[0])):
            organize_barentzwatch.organize_files()

        chunk = []
        for loc in location:

            chunk.append(self.load_barentswatch(self.path[data_label],loc))

        return xr.concat(chunk,'location')


    @staticmethod
    def load_barentswatch(path,location):

        filename  = '_'.join(['barentswatch',str(location)])+'.nc'

        open_data = xr.open_dataset(path+filename)

        return open_data

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
