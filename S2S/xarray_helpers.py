import xarray as xr
import numpy as np
import pandas as pd

def clim_mean(x):
    result = []
    for n in range(x.shape[0]):
        result.append(np.nanmean(np.delete(x,n,axis=0)))
    return np.array(result)

def clim_std(x):
    result = []
    for n in range(x.shape[0]):
        result.append(np.nanstd(np.delete(x,n,axis=0)))
    return np.array(result)

def clim_mean_c(x):
    return np.nanmean(x)

def clim_std_c(x):
    return np.nanstd(x)

def o_climatology(da,dim):
    print('\t performing xarray_helpers.o_climatology()')
    dim_name = dim.split('.')[0]
    mean,std = [],[]
    for label,data in list(da.groupby(dim)):
        mean.append(
                    xr.apply_ufunc(
                        clim_mean, data,
                        input_core_dims  = [[dim_name]],
                        output_core_dims = [[dim_name]],
                        vectorize=True, dask='parallelized'
                    )
                )
        std.append(
                    xr.apply_ufunc(
                        clim_std, data,
                        input_core_dims  = [[dim_name]],
                        output_core_dims = [[dim_name]],
                        vectorize=True, dask='parallelized'
                    )
                )

    return xr.concat(mean,dim_name).sortby(dim_name),\
                xr.concat(std,dim_name).sortby(dim_name)

def c_climatology(da,dim):
    print('\txarray_helpers.c_climatology()')
    dim_name = dim.split('.')[0]
    subgroup = dim.split('.')[1]
    mean,std = [],[]
    subgroup_label = []
    for label,data in list(da.groupby(dim)):
        subgroup_label.append(label)
        data = data.unstack()
        mean.append(
                    xr.apply_ufunc(
                        clim_mean_c, data,
                        input_core_dims  = [['time','member']],
                        output_core_dims = [[]],
                        vectorize=True, dask='parallelized'
                    )
                )
        std.append(
                    xr.apply_ufunc(
                        clim_std_c, data,
                        input_core_dims  = [['time','member']],
                        output_core_dims = [[]],
                        vectorize=True, dask='parallelized'
                    )
                )

    mean = xr.concat(mean,pd.Index(subgroup_label,name=subgroup))
    std  = xr.concat(std,pd.Index(subgroup_label,name=subgroup))
    return mean,std

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
    print('\txarray_helpers.climatology_to_validation_time()')
    time = validation_time.time
    step = validation_time.step

    tout = []
    for t in time:
        sout = []
        for s in step:
            sout.append(da.sel(month=pd.Timestamp((t+s).values).month,step=s))
        tout.append(xr.concat(sout,'step').drop('month'))
        p = 1-((time[-1]-t)/(time[-1]-time[0])).values
        print(str(round(p*100,1))+' % done',end='\r')
    return xr.concat(tout,time)

def assign_validation_time(ds):
    return ds.assign_coords(validation_time=ds.time+ds.step)

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

################################################################################
#######################  FUNCTIONS BELOW ARE DEPRICATED  #######################
################################################################################
def c_by_vt(da):
    out = []
    # da.expand_dims('validation_time')\
    #         .assign_coords(validation_time=da.time+da.step)
    for label,data in list(da.groupby('time.dayofyear')):
        out.append(
                    xr.apply_ufunc(
                        standard, data.unstack(),
                        input_core_dims  = [['time','member','step']],
                        output_core_dims = [['time','member','step']],
                        vectorize=True
                    ).rename('mean')
                )

    return xr.concat(out,'time').sortby(['time','step'])

def c_standardize_D(da,dim):

    print('\txarray_helpers.c_standardize()')
    dim_name = dim.split('.')[0]
    out = []
    for label,data in list(da.groupby(dim)):

        stack_name = list(data.dims)[-1]
        stack_dims = tuple(stack_name.split('_')[1:])
        stack_idx  = data[stack_name]

        data = data.unstack()
        out.append(
                    xr.apply_ufunc(
                        standard, data,
                        input_core_dims  = [['time','member']],
                        output_core_dims = [['time','member']],
                        vectorize=True
                    ).stack({stack_name:stack_dims})\
                        .assign_coords({stack_name:stack_idx})
                )
    return xr.concat(out,stack_name).unstack()

def c_standardize(da,dim):

    print('\txarray_helpers.c_standardize()')
    dim_name = dim.split('.')[0]
    out = []
    for label,data in list(da.groupby(dim)):

        data = data.unstack()
        out.append(
                    xr.apply_ufunc(
                        standard, data,
                        input_core_dims  = [['time','member']],
                        output_core_dims = [['time','member']]
                    )
                )
        print()
    return xr.concat(out,stack_name).unstack()

def ca_standardize(da):
    stack_name = list(da.dims)[-1]
    stack_dims = tuple(stack_name.split('_')[1:])
    da = da.unstack()

    out = xr.apply_ufunc(
        standard, da,
        input_core_dims  = [['time','member']],
        output_core_dims = [['time','member']],
        vectorize=True
    )
    return out.stack({stack_name:stack_dims})

def standard(x):
    return (x-np.nanmean(x))/np.nanstd(x)

def get_groups(ds,dim,keys):
    dim_name     = dim.split('.')[0]
    groups       = []
    grouped_data = ds.groupby(dim).groups
    for key in keys:
        idx = grouped_data[np.int64(key.values)]
        groups.append(ds.isel({dim_name:idx}))
    return xr.concat(groups,dim_name)

def keep_groups_of(ds,dim,members=5):
    dim_name     = dim.split('.')[0]
    groups       = []
    grouped_data = ds.groupby(dim).groups
    for key in list(grouped_data):
        idx = grouped_data[key]
        if len(idx) < members:
            pass
        else:
            groups.append(ds.isel({dim_name:idx}))
    return xr.concat(groups,dim_name).sortby(dim_name)

def match_times(cast,obs):
    """
    Match cast to obs and stack obs like cast

    args:
        cast: xarray.Dataset with time and step dimension
        obs: xarray.Dataset with time dimension

    returns
        obs: xarray.Dataset with time and step dimension
    """
    print('\t performing xarray_helpers.match_times()')

    step  = cast.step.sortby('step').to_pandas()
    obs   = obs.sortby('time')

    start = obs.time[0].to_pandas()
    end   = obs.time[-1].to_pandas()-step.max()

    cast = cast.sortby('time').sel(
                                time=slice(
                                    start,
                                    end
                                )
                            )
    out_obs  = []
    obs = obs.reindex(
            time=pd.date_range(
                start=start,
                end=end+step.max(),
                freq=step[2]-step[1]
            )
        )

    for time in cast.time:
        out_obs.append(obs.sel(time=time+cast.step).drop('time'))
    return cast,xr.concat(out_obs,cast.time)



def match_core(obs,times,steps):
    """
    args:
        obs: xarray.Dataset with time dimension
        times: xarray.DataArray
        steps: xarray.DataArray
    returns
        obs: xarray.Dataset with time and step dimension
    """
    print('\t\txarray_helpers.match_core()')

    obs     = obs.sortby('time')
    steps   = steps.sortby(steps)
    times   = times.sortby(times)
    times   = times.sel(time=slice(
                                obs.time.min(),
                                obs.time.max()-steps.max()
                                )
                            )
    out = []
    for time in times:
        out.append(obs.sel(time=time+steps).drop('time'))
    return xr.concat(out,times).sortby('time')

def stack_model(clim,cast):
    """
    Stack model like cast
    """
    print('\t performing xarray_helpers.stack_model(), an incredibly slow routine')
    clim   = clim.sortby('dayofyear')
    clim   = clim.reindex(dayofyear=np.arange(1,367,1,dtype='int'))
    t_out  = []
    for time in cast.time:
        date = pd.to_datetime(
                        xr.DataArray(
                            time.variable
                        ).to_pandas()
                    )
        t_out.append(clim.sel(dayofyear=date.day_of_year))
    return xr.concat(t_out,cast.time)

def stack_clim(clim,cast):
    """
    Stack climatology like cast
    """
    print('\t performing xarray_helpers.stack_clim(), an incredibly slow routine')
    clim   = clim.sortby('month')
    clim   = clim.reindex(month=np.arange(1,13,1,dtype='int'))
    t_out  = []
    for time in cast.time:
        s_out = []
        for step in cast.step:
            date = pd.to_datetime(
                            xr.DataArray(
                                time.variable+step.variable
                            ).to_pandas()
                        )
            s_out.append(clim.sel(month=date.month))
        t_out.append(xr.concat(s_out,cast.step))
    return xr.concat(t_out,cast.time)

def stack_like(obs,cast):
    """
    Match cast to obs and stack obs like cast

    args:
        cast: xarray.Dataset with time and step dimension
        obs: xarray.Dataset with time dimension

    returns
        obs: xarray.Dataset with time and step dimension

    Difference from match_time() indicated by *
    """
    print('\t performing xarray_helpers.stack_like()')
    obs = obs.sortby('time')
    cast = cast.sortby('time')

    step  = cast.step.to_pandas()
    time  = obs.time.to_pandas()

    start = time.min()
    end   = time.max()

    cast = cast.sel(time=slice(start,end-step.max()))

    obs = obs.reindex(
            time=pd.date_range(
                start=start,
                end=end,
                freq='D'
            )
        )

    out_obs  = []
    dims  = cast.dims['step']
    steps = cast.sortby('step').step
    for time in cast.time:
        out_obs.append(xr.concat([obs.sel(time=time).drop('time')]*dims,steps))
    return xr.concat(out_obs,cast.time)
