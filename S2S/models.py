import numpy as np
import pandas as pd
import xarray as xr
import statsmodels.api as sm
import scipy.stats as stats

import S2S.xarray_helpers as xh

def running_regression_CV(x,index,window=30):
    """
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
            # regression

            ############
            ymean.append(np.nanmean(filtered_pool))
            ystd.append(np.nanstd(filtered_pool))

        mean.append(np.array(ymean))
        std.append(np.array(ystd))

    return np.stack(mean,axis=-1),np.stack(std,axis=-1)

def regression1D(x,y):
    """
    Returns fitted coeffs by OLS for model
        y = intercept + slope*x + residuals

    Keeps the

    args:
        x: np.array dim: (n,)
        y: np.array dim: (n,)

    returns:
        intercept,slope: np.array,np.array dim: (n,),(n,)
    """
    res1 = []
    res2 = []
    for n in range(len(x)):
        X = np.delete(x,n)
        X = np.stack([np.ones_like(X),X],axis=1)
        res = sm.OLS(np.delete(y,n),X).fit()
        res1.append(res.params[0])
        res2.append(res.params[1])
    return np.array(res1),np.array(res2)

def regression2D(x1,x2,y):
    """
    Returns fitted coeffs by OLS for model
        y = intercept +slope1*x1 + slope2*x2 + residuals

    Keeps the

    args:
        x1: np.array dim: (n,)
        x2: np.array dim: (n,)
        y : np.array dim: (n,)

    returns:
        intercept,slope1,slope2: np.array,np.array dim: (n,),(n,)
    """
    res1 = []
    res2 = []
    res3 = []
    for n in range(y.shape[0]):
        X1 = np.delete(x1,n,axis=0).flatten()
        X2 = np.delete(x2,n,axis=0).flatten()
        X0 = np.ones_like(X1)
        X = np.stack(
                        [
                            X0,
                            X1,
                            X2
                        ],axis=1
                    )
        Y = np.delete(y,n,axis=0).flatten()
        res = sm.OLS(Y, X, missing='drop').fit()
        res1.append(res.params[0])
        res2.append(res.params[1])
        res3.append(res.params[2])
    return np.array(res1),np.array(res2),np.array(res3)

def combo(obs_t0,model_t,obs_t,dim='validation_time.month'):
    """
    Not very elegant and not particularly cheap,
    TODO: RE-DO

    args:
        obs_t0:  xarray.DataArray
        model_t: xarray.DataArray
        obs_t:   xarray.DataArray

    returns:


    """
    print('\tmodels.combo()')

    dim_name = dim.split('.')[0]
    subgroup = dim.split('.')[1]

    model    = model_t
    model_t  = model.mean('member',skipna=True)

    if obs_t0.dims != model.dims:
        obs_t0 = obs_t0.broadcast_like(model_t)

    if obs_t.dims != model.dims:
        obs_t  = obs_t.broadcast_like(model_t)

    ds = xr.merge(
                    [
                        obs_t0.rename('obs_t0'),
                        model_t.rename('model_t'),
                        obs_t.rename('obs_t')
                ],join='inner',compat='override'
            )

    groups = list(ds.groupby(dim))

    n = 0
    N = len(groups)

    xh.print_progress(n,N)
    e1,e2 = ' s','sl'

    stacked_dim = list(groups[0][1].dims)[-1]
    stack_dim_1 = stacked_dim.split('_')[-2]
    stack_dim_2 = stacked_dim.split('_')[-1]

    subgroup_data,subgroup_label = [],[]
    for n,(label,data) in enumerate(groups):

        data = data.unstack()

        tup = xr.apply_ufunc(
                regression2D, data.obs_t0, data.model_t, data.obs_t,
                input_core_dims  = [['time'],['time'],['time']],
                output_core_dims = [['time'],['time'],['time']],
                vectorize=True
                )

        subgroup_data.append(xr.merge(
            [
                tup[0].rename('intercept'),
                tup[1].rename('model_t'),
                tup[2].rename('obs_t0')
                ]
            ).stack({stacked_dim:(stack_dim_1,stack_dim_2)})
        )

        e1 += 'o'
        e2 += 'o'
        xh.print_progress(n+1,N,e=e1+' '+e2+'w')

    return xr.concat(subgroup_data,stacked_dim).unstack()


def clim_fc(mean,std,r=1,number_of_members=11):
    """
    Combines mean and std along a member dimension

    args:
        mean: xarray.DataArray
        std:  xarray.DataArray
    returns:
        clim_fc: xarray.DataArray with additional member dimension
    """
    mean = mean.expand_dims('member')
    std  = std.expand_dims('member')

    return xr.concat([mean-r*std,mean+r*std],'member')

def deterministic_gaussian_forecast(mean,std):
    """
    Generate random deterministic forecast of same dimensions as mean and std.
    Mean and std must have identical shape/dims.

    args:
        mean: xarray.DataArray
        std: xarray.DataArray

    returns:
        forecast: xarray.DataArray
    """

    return xr.apply_ufunc(
            np.random.normal, mean, std,
            input_core_dims  = [[],[]],
            output_core_dims = [[]],
            vectorize=True,dask='parallelized'
            )

def bias_adjustment_torrabla(forecast,observations,clim_std=None,window=30):
    """
    Forecsat calibration as of Eq. 2-4 in Torralba et al. (2017).

    References:

    Torralba, V., Doblas-Reyes, F. J., MacLeod, D., Christel, I., & Davis, M.
    (2017). Seasonal Climate Prediction:
    A New Source of Information for the Management of Wind Energy Resources,
    Journal of Applied Meteorology and Climatology, 56(5), 1231-1247.
    Retrieved Jun 22, 2021, from
    https://journals.ametsoc.org/view/journals/apme/56/5/jamc-d-16-0204.1.xml
    """

    x = forecast.mean('member')
    z = forecast - x

    if clim_std is None or window!=30:
        _,clim_std = xh.o_climatology(observations,window=window)

    ds = xr.merge(
                    [
                        forecast.rename('fc'),
                        observations.rename('obs')
                    ],join='inner',compat='override'
                )

    ds = xh.unstack_time(ds)

    rho = xr.apply_ufunc(
                correlation_CV,ds.fc.mean('member'),ds.obs,ds.dayofyear,window,
                input_core_dims  = [
                                    ['year','dayofyear'],
                                    ['year','dayofyear'],
                                    ['dayofyear'],
                                    []
                                ],
                output_core_dims = [['year','dayofyear']],
                vectorize=True
    )

    rho = xh.stack_time(rho)

    sigma_ens = xr.apply_ufunc(
                std,ds.fc.mean('member'),ds.dayofyear,window,
                input_core_dims  = [
                                    ['year','dayofyear'],
                                    ['dayofyear'],
                                    []
                                ],
                output_core_dims = [['year','dayofyear']],
                vectorize=True
    )
    sigma_ens = xh.stack_time(sigma_ens)

    sigma_ref = clim_std

    sigma_e   = fc.std('member')


    alpha = xr.ufuncs.fabs(rho) * sigma_ref/sigma_ens

    beta  = ( 1 - )

def correlation_CV(x,y,index,window=30):
    """
    Computes correlation of x against y (rho), keeping dim -1 and -2.
    Dim -1 must be 'dayofyear', with the corresponding days given in index.
    Dim -2 must be 'year'.

    args:
        x:      np.array of float, with day of year as index -1 and year as
                index -2
        index:  np.array of int, 1-dimensional holding dayofyear corresponding
                to dim -1 of x

    returns
        rho:   np.array of float, with day of year as index -1 and year as
                dim -2

    dimensions requirements:
        name            dim

        year            -2
        dayofyear       -1
    """
    rho   = []

    pad   = window//2

    x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]
    y     = np.pad(y,pad,mode='wrap')[pad:-pad,:]

    index = np.pad(index,pad,mode='wrap')

    index[-pad:] += index[-pad-1]
    index[:pad]  -= index[-pad-1]

    for ii,idx in enumerate(index[pad:-pad]):

        # pool all values that falls within window
        xpool = x[...,np.abs(index-idx)<=pad]
        ypool = y[...,np.abs(index-idx)<=pad]

        yrho = []
        for yy in range(xpool.shape[-2]):

            # delete the relevant year from pool (for cross validation)
            filtered_xpool = np.delete(xpool,yy,axis=-2).flatten()
            filtered_ypool = np.delete(ypool,yy,axis=-2).flatten()

            idx_bool = ~np.logical_or(
                                np.isnan(filtered_xpool),
                                np.isnan(filtered_ypool)
                            )

            r,p = stats.pearsonr(
                                    filtered_xpool[idx_bool],
                                    filtered_ypool[idx_bool]
                                )

            yrho.append(r)

        rho.append(np.array(yrho))

    return np.stack(rho,axis=-1)

def std(x,index,window=30):
    """
    Computes std of x, keeping dim -1 and -2.
    Dim -1 must be 'dayofyear', with the corresponding days given in index.
    Dim -2 must be 'year'.

    args:
        x:      np.array of float, with day of year as index -1 and year as
                index -2
        index:  np.array of int, 1-dimensional holding dayofyear corresponding
                to dim -1 of x

    returns
        std:   np.array of float, with day of year as index -1 and year as
                dim -2

    dimensions requirements:
        name            dim

        year            -2
        dayofyear       -1
    """
    std   = []

    pad   = window//2

    x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]
    index = np.pad(index,pad,mode='wrap')

    index[-pad:] += index[-pad-1]
    index[:pad]  -= index[-pad-1]

    for ii,idx in enumerate(index[pad:-pad]):

        # pool all values that falls within window
        pool = x[...,np.abs(index-idx)<=pad]

        std.append(np.full_like(pool[...,0],np.nanstd(pool)))

    return np.stack(std,axis=-1)

################################################################################
############################# Depricated #######################################
################################################################################
def persistence(predictor,response,var,dim='time.dayofyear'):
    """
    Not very elegant and probably not particularly cheap
    """
    dim_name = dim.split('.')[0]
    print('\t performing models.persistence()')
    ds = xr.merge(
                    [
                        predictor.rename({var:'predictor'}),
                        response.rename({var:'response'})
                ],join='inner',compat='override'
            )
    n = 0
    N = len(list(ds.groupby(dim)))
    data_out = []
    for label,data in list(ds.groupby(dim)):

        tup = xr.apply_ufunc(
                regression1D, data.predictor, data.response,
                input_core_dims  = [[dim_name], [dim_name]],
                output_core_dims = [[dim_name], [dim_name]],
                vectorize=True
                )
        data_out.append(xr.merge(
            [
                tup[0].rename('intercept'),
                tup[1].rename('slope')
                ]
            )
        )

    return xr.concat(data_out,dim_name).sortby(dim_name)

def persistence2(predictor,response,var):
    """
    Not very elegant and probably not particularly cheap
    """
    print('\t performing models.persistence2()')
    ds = xr.merge(
                    [
                        predictor.rename({var:'predictor'}),
                        response.rename({var:'response'})
                ],join='inner',compat='override'
            )
    n = 0
    N = len(list(ds.groupby('time.dayofyear')))
    date_out,date_label_out = [],[]
    for date_label,date_group in list(ds.groupby('time.dayofyear')):
        n += 1
        print('\t\t group ',n,' of ',N)
        time_out = []
        for time in date_group.time:
            data = date_group.sel(
                            time=date_group.time\
                                .where(date_group.time!=time,drop=True
                                )
                            )
            time_out.append(xr.apply_ufunc(
                    regression1D2, data.predictor, data.response,
                    input_core_dims  = [['time'], ['time']],
                    output_core_dims = [[]],
                    vectorize=True, dask='parallelized'
                    ).rename('slope')
                )
        date_out.append(xr.concat(time_out,date_group.time))
    return xr.concat(date_out,'time')

def combo2(observation,model,response,var):
    """
    Not very elegant and probably not particularly cheap
    """
    print('\t performing models.combo2(), takes a while')
    observation = xr.concat([observation]*model.dims['member'],model.member)
    response    = xr.concat([response]*model.dims['member'],model.member)

    ds = xr.merge(
                    [
                        observation.rename({var:'observation'}),
                        model.rename({var:'model'}),
                        response.rename({var:'response'})
                ],join='inner',compat='override'
            )
    N = len(list(ds.groupby('time.dayofyear')))
    n = 0
    date_out,date_label_out = [],[]
    for date_label,date_group in list(ds.groupby('time.dayofyear')):
        time_out = []
        n += 1
        print('\t\t group ',n,' of ',N)
        for time in date_group.time:
            data = date_group.sel(
                            time=date_group.time\
                                .where(date_group.time!=time,drop=True
                                )
                            )
            tup = xr.apply_ufunc(
                    regression2D2, data.observation, data.model, data.response,
                    input_core_dims  = [['time'], ['time'], ['time']],
                    output_core_dims = [[],[]],
                    vectorize=True
                    )

            time_out.append(xr.merge(
                [
                    tup[0].rename('slope_obs'),
                    tup[1].rename('slope_model')
                    ]
                )
            )
        date_out.append(xr.concat(time_out,date_group.time))
    return xr.concat(date_out,'time')

def regression1D2(x,y):
    result = sm.OLS(y, x).fit()
    return result.params[0]

def regression2D2(x1,x2,y):
    X = np.stack([x1,x2],axis=1)
    result = sm.OLS(y, X).fit()
    return result.params[0],result.params[1]

def reg_m(y, x):
    """
    Same as Erik's eide.py - ewkutilis.py, but with reversed order of predictors
    and therefore also output params
    """
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((ones,x[0])))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((X,ele)))
    results = sm.OLS(y, X).fit()
    return results

def a_and_b(model):
    return xr.merge(
                [
                    xr.DataArray(
                        data = model.params[0],name='intercept'
                        ),
                    xr.DataArray(
                        data = model.params[1],name='slope'
                    )
                ]
            )

def a_b_and_c(model):
    return xr.merge(
                [
                    xr.DataArray(
                        data = model.params[0],name='intercept'
                        ),
                    xr.DataArray(
                        data = model.params[1],name='slope_obs'
                    ),
                    xr.DataArray(
                        data = model.params[2],name='slope_model'
                    )
                ]
            )
