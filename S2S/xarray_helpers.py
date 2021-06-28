import xarray as xr
import numpy as np
import pandas as pd
import scipy.stats as stats

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

def running_clim(x,index,window=30):
    """
    Computes mean and standard deviation over x keeping dim -1. Dim -1 must be
    'dayofyear', with the corresponding days given in index.

    args:
        x:      np.array of float, with day of year as index -1
        index:  np.array of int, 1-dimensional holding dayofyear corresponding
                to dim -1 of x

    returns
        mean:   np.array of float, with day of year as index -1 and year as
                dim -2
        std:   np.array of float, with day of year as index -1 and year as
                dim -2

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

        # pool all values that falls within window
        pool = x[...,np.abs(index-idx)<=pad]

        if np.isfinite(pool).sum() > 1:
            mean.append(np.full_like(pool[0][...,0],np.nanmean(pool)))
            std.append(np.full_like(pool[0][...,0],np.nanstd(pool)))
        else:
            mean.append(np.full_like(pool[0][...,0],np.nan))
            std.append(np.full_like(pool[0][...,0],np.nan))

    return np.stack(mean,axis=-1),np.stack(std,axis=-1)

def running_clim_CV(x,index,window=30):
    """
    Like running_clim() but leaves out current year in computation for
    cross validation

    Computes mean and standard deviation over x keeping dim -1 and -2.
    Dim -1 must be 'dayofyear', with the corresponding days given in index.
    Dim -2 must be 'year'.

    args:
        x:      np.array of float, with day of year as index -1 and year as
                index -2
        index:  np.array of int, 1-dimensional holding dayofyear corresponding
                to dim -1 of x

    returns
        mean:   np.array of float, with day of year as index -1 and year as
                dim -2
        std:   np.array of float, with day of year as index -1 and year as
                dim -2

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

        # pool all values that falls within window
        pool = x[...,np.abs(index-idx)<=pad]

        ymean,ystd = [],[]
        for yy in range(pool.shape[-2]):

            # delete the relevant year from pool (for cross validation)
            filtered_pool = np.delete(pool,yy,axis=-2)

            if np.isfinite(filtered_pool).sum() > 1:
                ymean.append(np.nanmean(filtered_pool))
                ystd.append(np.nanstd(filtered_pool))
            else:
                ymean.append(np.nan)
                ystd.append(np.nan)

        mean.append(np.array(ymean))
        std.append(np.array(ystd))

    return np.stack(mean,axis=-1),np.stack(std,axis=-1)

def o_climatology(da,window=30):
    """
    Climatology with centered initialization time, using cross validation
    using a 30-day window.

    args:
        da: xarray.DataArray with 'time' dimension

    returns:
        mean: xarray.DataArray, like da
        std: xarray.DataArray, like da

    time: datetime-like
    """

    print('\t\txarray_helpers.o_climatology()')

    da = unstack_time(da)

    # to all year,dayofyear matrices in da, apply runnning_clim_CV
    mean,std = xr.apply_ufunc(
            running_clim_CV, da, da.dayofyear,window,
            input_core_dims  = [['year','dayofyear'],['dayofyear'],[]],
            output_core_dims = [['year','dayofyear'],['year','dayofyear']],
            vectorize=True
        )

    # re-assing time dimension to da from year,dayofyear
    return stack_time(mean),stack_time(std)

def c_climatology(da):
    """
    Climatology with centered initialization time
    """

    print('\t\txarray_helpers.c_climatology()')

    da = unstack_time(da)

    # to all year,dayofyear matrices in da, apply runnning_clim
    mean,std = xr.apply_ufunc(
            running_clim, da, da.dayofyear,
            input_core_dims  = [['member','year','dayofyear'],['dayofyear']],
            output_core_dims = [['year','dayofyear'],['year','dayofyear']],
            vectorize=True
        )

    # re-assing time dimension to da from year,dayofyear
    return stack_time(mean),stack_time(std)

def unstack_time(da):
    """
    Splits time dimension in da into a 'year' and a 'dayofyear' dimension.
    Coordinate time must be datetime-like.

    args:
        da: xarray.DataArray, requires dims: time

    returns:
        da: xarray.DataArray, new dimensions: year, dayofyear
    """
    da   = da.sortby('time')
    time = da.time

    # create an mulitindex mapping time -> (year,dayofyear)
    stacked_time = pd.MultiIndex.from_arrays(
                                            [
                                                da.time.dt.year.to_pandas(),
                                                da.time.dt.dayofyear.to_pandas()
                                            ],names=('year','dayofyear')
                                        )

    # re-assing time from datetime like to multiindex (year,dayofyear) and
    # and split mutliindex into year and dauofyear dimension
    return da.assign_coords(time=stacked_time).unstack('time')

def stack_time(da):
    """
    Stacks 'year' and 'dayofyear' dimensions in a xarray.DataArray to a 'time'
    dimension.

    args:
        da: xarray.DataArray, requires dims: year and dayofyear

    returns:
        da: xarray.DataArray, new dimension: time
    """

    da = da.stack(time=('year','dayofyear'))

    time = []
    for year,dayofyear in zip(da.year.values,da.dayofyear.values):
        year       = pd.to_datetime(year, format='%Y')
        dayofyear  = pd.Timedelta(dayofyear-1,'D')
        time.append(year+dayofyear)

    da = da.assign_coords(time=time)

    return da

def assign_validation_time(ds):
    """
    Add validation_time coordinates to xarray.Dataset/DataArray with 'time' and
    'step' dimensions.

    validation_time = time + step

    time: datetime-like
    step: pd.Timedelta

    args:
        ds: xarray.Dataset/DataArray with 'time' and 'step' dimensions

    returns:
        ds: xarray.Dataset/DataArray with 'validation_time' dimension
    """
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

        try:
            out.append(
                xr.open_dataset(config['VALID_DB']+'/'+loc+'_'+filename+'.nc')
                )
        except FileNotFoundError:
            pass

    return xr.concat(out,'location')

def xtrapolate_NAN(x):
    """
    Underfunction of extrapolate_land_mask

    Inter/extrapolate NaN values by meaning over the sorrounding grid points
    ignoring NaN values.

    if x at grid point [i,j] equals NaN, if x[i,j] == NaN, then

        x[i,j] = np.nanmean( x[i-1,j], x[i+1,j], x[i,j-1], x[i,j+1] )

    if all x[i-1,j], x[i+1,j], x[i,j-1], x[i,j+1] are NaN values, then
    x[i,j] = NaN

    args:
        x: np.array 2-dimensional

    returns:
        x: np.array 2-dimensional

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

# def isolate_highest_r(x,y,index,window=30):
#     """
#     """
#     # Flatten grid
#     y = y.reshape(y.shape[0],y.shape[1],y.shape[2],-1)
#
#     # keep y with member dim
#     Y = y
#
#     # mean over member dim
#     y = y.mean(2)
#
#     pad = window//2
#
#     x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]
#     y     = np.pad(y,pad,mode='wrap')[pad:-pad,:,pad:-pad]
#
#     index = np.pad(index,pad,mode='wrap')
#
#     index[-pad:] += index[-pad-1]
#     index[:pad]  -= index[-pad-1]
#
#     ys = [] # for each dayofyear
#     for ii,idx in enumerate(index[pad:-pad]):
#
#         # pool all values that falls within window around respective dayofyear
#         xpool = x[:,np.abs(index-idx)<=pad]
#         ypool = y[:,np.abs(index-idx)<=pad]
#
#         # flatten year and dayofyear
#         xpool = xpool.reshape(-1)
#         ypool = ypool.reshape(-1,ypool.shape[-1])
#
#         r = [] # for each gridpoint
#         for ii_y in range(ypool.shape[-1]):
#
#             # keep only nonnan forecast-observation pairs
#             idx_bool = ~np.logical_or(
#                                 np.isnan(xpool),
#                                 np.isnan(ypool[:,ii_y])
#                             )
#
#             xp = xpool[idx_bool]
#             yp = ypool[idx_bool,ii_y]
#
#             # if all nan (most likely landmask) give -99 as correlation coef.
#             if idx_bool.sum()==0:
#                 r.append(-99)
#             # otherwise pearsons r
#             else:
#                 r.append(stats.pearsonr(xp,yp)[0])
#
#         r = np.array(r)
#
#         # pick the gridpoint of highest correlation to represent the observation
#         # at the respective dayofyear
#         ys.append(Y[:,ii,:,np.argmax(r)])
#
#     # stack along dayofyear dimension
#     return np.stack(ys,axis=1)

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
    """
    Apply runnning mean on hindcast per location. Labels centered in the window.

    args:
        hindcast: xarray.DataArray with step dimension
        window:   number of steps included in window

    returns:
        hindcast: xarray.DataArray
    """

    N   = len(hindcast.location)
    out = []

    for n,loc in enumerate(hindcast.location):

        print_progress(n,N)

        h = hindcast.sel(location=loc)

        out.append(h.rolling(step=window,center=True)\
                .mean().dropna('step'))
    return xr.concat(out,'location')

def absolute(u,v):
    """
    Compute the absoulte value of two horizontal components.

    U = sqrt( u**2 + v**2 )

    args:
        u: np.array n-dimensional
        v: np.array n-dimensional

    returns:
        U: np.array n-dimensional
    """
    return np.sqrt( u**2 + v**2 )

def cor_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    Written by stackoverflow user: dbliss,
    Dowloaded June 25 2021 from https://stackoverflow.com/a/30145770
    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                  mu_y[np.newaxis, :])
    return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
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
