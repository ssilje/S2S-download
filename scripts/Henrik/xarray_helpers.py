import xarray as xr
import numpy as np
import pandas as pd

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
    return xr.concat(groups,dim_name)

def match_times(cast,obs):
    """
    Match cast to obs and stack obs like cast

    args:
        cast: xarray.Dataset with time and step dimension
        obs: xarray.Dataset with time dimension

    returns
        obs: xarray.Dataset with time and step dimension
    """
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
                end=end,
                freq=step[2]-step[1]
            )
        )

    for time in cast.time:
        out_obs.append(obs.sel(time=time+cast.step).drop('time'))
    return cast,xr.concat(out_obs,cast.time)

def stack_climatology(clim,cast):
    """
    Stack climatology like cast

    args:
        cast: xarray.Dataset with time and step dimension
        clim: xarray.Dataset with day/month-ofyear dimension

    returns
        clim: xarray.Dataset with time and step dimension
    """
    clim   = clim.sortby('dayofyear')
    clim   = clim.reindex(dayofyear=np.arange(1,367,1,dtype='int'))
    t_out  = []
    for time in cast.time:
        s_out = []
        for step in cast.sel(time=time).step:
            date = pd.to_datetime(
                        xr.DataArray(
                            time.variable+step.variable
                        ).to_pandas()
                    )

            s_out.append(clim.sel(dayofyear=date.day_of_year))
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
    cast  = cast.sortby('step')
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
                end=end,
                freq=step[2]-step[1]
            )
        )

    steps = cast.sortby('step').step
    for time in cast.time:
        step_temp = []
        for s in steps:
            step_temp.append(obs.sel(time=time).drop('time'))
        out_obs.append(xr.concat(step_temp,steps)) #*
    return xr.concat(out_obs,cast.time)
