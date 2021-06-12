import xarray as xr
import numpy as np
import pandas as pd

from S2S.local_configuration import config

from S2S.graphics import graphics as gr

def print_progress(n,N,i='',e=''):
    """
    Print progression

    args:
        n: float or int, current step
        N: float or int, total number of steps
    """
    print('\t\t'+i+str(round((n/N)*100,1))+' % done'+e,end='\r')

def running_clim_CV(x,index,window=30):
    """
    Like runing_clim() but leaves out current year in computation for
    cross validation

    dimensions requirements:
        name            dim

        year            -2
        dayofyear       -1
    """
    mean  = []
    std   = []

    pad   = window//2

    x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]
    index = np.pad(index,pad,mode='wrap')

    index[-pad:] += index[-pad-1]
    index[:pad]  -= index[-pad-1]

    for ii,idx in enumerate(index[pad:-pad]):

        pool = x[...,np.abs(index-idx)<=pad]

        ymean,ystd = [],[]
        for yy in range(pool.shape[-2]):

            filtered_pool = np.delete(pool,yy,axis=-2)
            ymean.append(np.nanmean(filtered_pool))
            ystd.append(np.nanstd(filtered_pool))

        mean.append(np.array(ymean))
        std.append(np.array(ystd))

    return np.stack(mean,axis=-1),np.stack(std,axis=-1)

def o_climatology(da):
    """
    Climatology with centered initialization time, using cross validation
    """

    print('\t\txarray_helpers.o_climatology()')

    da   = da.sortby('time')
    time = da.time
    stacked_time = pd.MultiIndex.from_arrays(
                                            [
                                                da.time.dt.year.to_pandas(),
                                                da.time.dt.dayofyear.to_pandas()
                                            ],names=('year','dayofyear')
                                        )

    da = da.assign_coords(time=stacked_time).unstack('time')

    mean,std = xr.apply_ufunc(
            running_clim_CV, da, da.dayofyear,
            input_core_dims  = [['year','dayofyear'],['dayofyear']],
            output_core_dims = [['year','dayofyear'],['year','dayofyear']],
            vectorize=True
        )

    return stack_time(mean),stack_time(std)

def running_clim(x,index,window=30):
    """
    dimensions requirements:
        name            dim

        dayofyear      -1

        member          0
    """
    mean  = []
    std   = []

    pad   = window//2

    x     = np.pad(x,pad,mode='wrap')[pad:-pad,pad:-pad,:]
    index = np.pad(index,pad,mode='wrap')

    index[-pad:] += index[-pad-1]
    index[:pad]  -= index[-pad-1]

    for ii,idx in enumerate(index[pad:-pad]):

        pool = x[...,np.abs(index-idx)<=pad]

        mean.append(np.full_like(pool[0][...,0],np.nanmean(pool)))
        std.append(np.full_like(pool[0][...,0],np.nanstd(pool)))

    return np.stack(mean,axis=-1),np.stack(std,axis=-1)

def c_climatology(da):
    """
    Climatology with centered initialization time
    """

    print('\t\txarray_helpers.c_climatology()')

    da   = da.sortby('time')
    time = da.time
    stacked_time = pd.MultiIndex.from_arrays(
                                            [
                                                da.time.dt.year.to_pandas(),
                                                da.time.dt.dayofyear.to_pandas()
                                            ],names=('year','dayofyear')
                                        )

    da = da.assign_coords(time=stacked_time).unstack('time')

    mean,std = xr.apply_ufunc(
            running_clim, da, da.dayofyear,
            input_core_dims  = [['member','year','dayofyear'],['dayofyear']],
            output_core_dims = [['year','dayofyear'],['year','dayofyear']],
            vectorize=True
        )

    return stack_time(mean),stack_time(std)

def stack_time(da):

    da = da.stack(time=('year','dayofyear'))

    time = []
    for year,dayofyear in zip(da.year.values,da.dayofyear.values):
        year       = pd.to_datetime(year, format='%Y')
        dayofyear  = pd.Timedelta(dayofyear-1,'D')
        time.append(year+dayofyear)

    da = da.assign_coords(time=time)

    return da

def assign_validation_time(ds):
    return ds.assign_coords(validation_time=ds.time+ds.step)

def store_by_location(da,filename):
    """
    Store file with location dim, per location

    args:
        da: xarray.DataArray or xarray.Dataset, dims required: location
        filename: str
    """
    print('\txarray_helpers.store_by_location()')

    if not hasattr(da.location,'__iter__'):
        loc = da.location
        da.to_netcdf(config['VALID_DB']+'/'+loc+'_'+filename+'.nc')

    else:
        N   = len(da.location)
        out = []

        for n,loc in enumerate(da.location.values):

            print_progress(n,N)

            da.sel(location=loc)\
                .to_netcdf(config['VALID_DB']+'/'+loc+'_'+filename+'.nc')

def load_by_location(location,filename):
    """
    Load file with location in in filename and concatinate along location dim

    args:
        location: str
        filename: str

    returns:
        da: xarray.DataArray or xarray.Dataset
    """
    print('\txarray_helpers.load_by_location()')
    N   = len(location)
    out = []

    for n,loc in enumerate(location.values):

        print_progress(n,N)

        out.append(
            xr.open_dataset(config['VALID_DB']+'/'+loc+'_'+filename+'.nc')
            )

    return xr.concat(out,'location')

def xtrapolate_NAN(x):
    """
    Underfunction of extrapolate_land_mask

    TODO: comment
    """

    x = np.pad(x,pad_width=1,mode='constant',constant_values=np.nan)

    x_int = np.nanmean(
                    np.stack(
                        [x[1:-1,:-2],x[1:-1,2:],x[:-2,1:-1],x[2:,1:-1]]
                ),axis=0
            )

    idx = np.isnan(x)[1:-1,1:-1]

    x.setflags(write=1)

    x[1:-1,1:-1][idx] = x_int[idx]

    x = x[1:-1,1:-1]

    # continue until grid is full
    if np.isnan(np.sum(x)):
        x = xtrapolate_NAN(x)

    return x

def extrapolate_land_mask(da):
    """
    Fill/extrapolate NaN values by linearly interpolating neighbouring values.
    Routine is repeated until all NaN (except boundaries) are filled in.

    args:
        da: xarray.DataArray or xarray.Dataset, dims required: lat,lon

    returns:
        da: NaNs filled in by linearly interpolating neighbours
    """
    return xr.apply_ufunc(
        xtrapolate_NAN, da,
        input_core_dims  = [['lon','lat']],
        output_core_dims = [['lon','lat']],
        vectorize=True
    )

def interp_to_loc(observations,hindcast):
    """
    Interpolate gridded hindcast to point locations specified by location
    dimension in observations. Returns hindcast with location dimension.

    args:
        observations: xarray.DataArray or xarray.Dataset, dims required: location,lat,lon
        hindcast: xarray.DataArray or xarray.Dataset, dims required: lat,lon

    returns:
        hindcast: xarray.DataArray or xarray.Dataset
    """
    N        = len(observations.location)
    out      = []

    for n,loc in enumerate(observations.location):

        print_progress(n,N)

        o = observations.sel(location=loc)

        lon = xr.ufuncs.fabs(hindcast.lon - o.lon)
        lat = xr.ufuncs.fabs(hindcast.lat - o.lat)

        h = hindcast.where(lon<2,drop=True)
        h = h.where(lat<2,drop=True)

        out.append(
                h.interp(
                    lon=o.lon,
                    lat=o.lat,
                    method='linear'
                    )
                )
    return xr.concat(out,'location')

################################################################################
################################################################################
################################################################################
def climatology_to_validation_time(da,validation_time):
    """
    TODO: Vectorize, incredibly slow routine but haven't figured out how to do
    this in a smart way yet...

    args:
        da:              xr.DataArray with step and month dimension
        validation_time: xr.DataArray with step + time dimension

    returns:
        climatology stacked to match validation time: xarray.DataArray
    """
    print('\t\txarray_helpers.climatology_to_validation_time()')
    time = validation_time.time
    step = validation_time.step

    tout = []
    for t in time:
        sout = []
        for s in step:
            # can switch pd.Timestamp((t+s).values).month with (t+s).dt.month (?)
            sout.append(da.sel(month=pd.Timestamp((t+s).values).month,step=s))
        tout.append(xr.concat(sout,'step').drop('month'))

        print_progress((t-(time[0])).values,(time[-1]-time[0]).values)

    return xr.concat(tout,time)



def at_validation(obs,vt,ddays=1):
    """
    args:
        obs:   xarray.Dataset with time dimension
        vt:    xarray.DataArray time + step dimension
        ddays: int, tolerance, in days, for matching times
    returns
        observations: xarray.Dataset with time and step dimension
    """
    print('\t\txarray_helpers.at_validation()')

    obs  = obs.sortby('time')
    vt   = vt.sortby(['time','step'])
    vt   = vt.sel(time=slice(obs.time.min(),obs.time.max()))
    time = vt.time
    step = vt.step

    out  = []
    for t in time:
        o = obs.reindex(
                indexers   = {'time':(t+step).data},
                method     = 'nearest',
                tolerance  = pd.Timedelta(ddays,'D'),
                fill_value = np.nan
                )
        o = o.assign_coords(time=o.time-t).rename(time='step')
        out.append(o)

    out = xr.concat(out,time).sortby(['time','step'])
    return assign_validation_time(out)



def running_mean(hindcast,window):

    N   = len(hindcast.location)
    out = []

    for n,loc in enumerate(hindcast.location):

        print_progress(n,N)

        h = hindcast.sel(location=loc)

        out.append(h.rolling(step=window,center=True)\
                .mean().dropna('step'))
    return xr.concat(out,'location')

def running_mean_o(observations,window):

    N   = len(observations.location)
    out = []

    for n,loc in enumerate(observations.location):

        print_progress(n,N)

        h = observations.sel(location=loc)

        out.append(h.rolling(time=window,center=True,min_periods=1)\
                .mean().dropna('time'))

    return xr.concat(out,'location',join='outer')

# def depricated_clim_c():
#     dim_name = dim.split('.')[0]
#     subgroup = dim.split('.')[1]
#     mean,std = [],[]
#     subgroup_label = []
#     for label,data in list(da.groupby(dim)):
#         subgroup_label.append(label)
#         data = data.unstack()
#         mean.append(
#                     xr.apply_ufunc(
#                         clim_mean_c, data,
#                         input_core_dims  = [['time','member']],
#                         output_core_dims = [[]],
#                         vectorize=True, dask='parallelized'
#                     )
#                 )
#         std.append(
#                     xr.apply_ufunc(
#                         clim_std_c, data,
#                         input_core_dims  = [['time','member']],
#                         output_core_dims = [[]],
#                         vectorize=True, dask='parallelized'
#                     )
#                 )
#
#     mean = xr.concat(mean,pd.Index(subgroup_label,name=subgroup))
#     std  = xr.concat(std,pd.Index(subgroup_label,name=subgroup))
#     return mean,std

# def depricated_o_climatology(da,dim):
#     print('\t\txarray_helpers.o_climatology()')
#     dim_name = dim.split('.')[0]
#     mean,std = [],[]
#     for label,data in list(da.groupby(dim)):
#         mean.append(
#                     xr.apply_ufunc(
#                         clim_mean, data,
#                         input_core_dims  = [[dim_name]],
#                         output_core_dims = [[dim_name]],
#                         vectorize=True, dask='parallelized'
#                     )
#                 )
#         std.append(
#                     xr.apply_ufunc(
#                         clim_std, data,
#                         input_core_dims  = [[dim_name]],
#                         output_core_dims = [[dim_name]],
#                         vectorize=True, dask='parallelized'
#                     )
#                 )
#
#     return xr.concat(mean,dim_name).sortby(dim_name),\
#                 xr.concat(std,dim_name).sortby(dim_name)
