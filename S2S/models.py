import numpy as np
import pandas as pd
import xarray as xr
import statsmodels.api as sm
import scipy.stats as stats

import S2S.xarray_helpers as xh

################################################################################
############################# xarray routines ##################################
################################################################################

def clim_fc(mean,std,r=1):
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

    return xr.concat([mean-r*std,mean+r*std],pd.Index([1,2],name='member'))

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

def bias_adjustment_torralba(
                                forecast,
                                observations,
                                clim_std=None,
                                window=30,
                                spread_only=False,
                                cluster_name=None
                            ):
    """
    Forecast calibration as of Eq. 2-4 in Torralba et al. (2017).

    References:

    Torralba, V., Doblas-Reyes, F. J., MacLeod, D., Christel, I., & Davis, M.
    (2017). Seasonal Climate Prediction:
    A New Source of Information for the Management of Wind Energy Resources,
    Journal of Applied Meteorology and Climatology, 56(5), 1231-1247.
    Retrieved Jun 22, 2021, from
    https://journals.ametsoc.org/view/journals/apme/56/5/jamc-d-16-0204.1.xml
    """

    if cluster_name:
        icd = [cluster_name,'year','dayofyear']
    else:
        icd = ['year','dayofyear']

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
                                    icd,
                                    icd,
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
                                    icd,
                                    ['dayofyear'],
                                    []
                                ],
                output_core_dims = [icd],
                vectorize=True
    )
    sigma_ens = xh.stack_time(sigma_ens)

    sigma_ref = clim_std

    sigma_e   = forecast.std('member')

    if spread_only:
        alpha = 1.
    else:
        alpha = xr.ufuncs.fabs(rho) * ( sigma_ref/sigma_ens )

    beta  = xr.ufuncs.sqrt( 1 - rho**2 ) * ( sigma_ref/sigma_e )

    y = alpha * x + beta * z

    return y

def persistence(init_value,observations,window=30):
    """
    Input must be anomlies.
    """

    print('\tmodels.persistence()')
    ds = xr.merge(
                    [
                        init_value.rename('iv'),
                        observations.rename('o')
                ],join='inner',compat='override'
            )

    ds  = xh.unstack_time(ds)

    rho = xr.apply_ufunc(
                correlation_CV,ds.iv,ds.o,ds.dayofyear,window,
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

    try:
        rho = rho.drop('validation_time')
    except (AttributeError, ValueError):
        pass

    return rho * init_value

def combo(
            init_value,
            model,
            observations,
            window=30,
            lim=1,
            sub=np.nan,
            cluster_name=None
        ):
    """
    Input must be anomalies.
    """

    if cluster_name:
        icd = [cluster_name,'year','dayofyear']
    else:
        icd = ['year','dayofyear']

    print('\t performing models.persistence()')
    ds = xr.merge(
                    [
                        init_value.rename('iv'),
                        model.rename('mod').mean('member'),
                        observations.rename('o')
                ],join='inner',compat='override'
            )

    ds  = xh.unstack_time(ds)

    alpha,beta = xr.apply_ufunc(
                running_regression_CV,
                ds.iv,
                ds.mod,
                ds.o,
                ds.dayofyear,
                window,
                lim,
                sub,
                input_core_dims  = [
                                    icd,
                                    icd,
                                    icd,
                                    ['dayofyear'],
                                    [],
                                    [],
                                    []
                                ],
                output_core_dims = [['year','dayofyear'],['year','dayofyear']],
                vectorize=True
            )

    try:
        alpha = alpha.drop('validation_time')
    except (AttributeError, ValueError):
        pass
    try:
        beta = beta.drop('validation_time')
    except (AttributeError, ValueError):
        pass

    alpha = xh.stack_time(alpha)
    beta  = xh.stack_time(beta)

    return alpha * init_value + beta * model.mean('member')

################################################################################
########################### NumPy routines #####################################
################################################################################
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

    Returns 2 dimensional array (year,dayofyear)
    """
    rho   = []

    pad   = window//2

    if len(x.shape)==2:
        x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]
        y     = np.pad(y,pad,mode='wrap')[pad:-pad,:]

    if len(x.shape)==3:
        x     = np.pad(x,pad,mode='wrap')[pad:-pad,pad:-pad,:]
        y     = np.pad(y,pad,mode='wrap')[pad:-pad,pad:-pad,:]

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

            idx_bool = np.logical_and(
                                np.isfinite(filtered_xpool),
                                np.isfinite(filtered_ypool)
                            )
            if idx_bool.sum()<2:
                r = np.nan

            else:
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
clim_fc = models.clim_fc(point_observations.mean,point_observations.std)
pers    = models.persistence(
                init_value   = point_observations.init_a,
                observations = point_observations.data_a
                )

graphics.timeseries(
                        observations    = point_observations.data_a,
                        cast            = [pers,point_hindcast.data_a],
                        lead_time       = [9,16],
                        clabs           = ['persistence','EC'],
                        filename        = 'BW_persistence',
                        title           = 'Barentswatch EC'
                    )

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

    Returns n dimensional array (...n-2,year,dayofyear) n is number of input dim
    """
    std   = []

    pad   = window//2

    if len(x.shape)==2:
        x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]

    if len(x.shape)==3:
        x     = np.pad(x,pad,mode='wrap')[pad:-pad,pad:-pad,:]

    index = np.pad(index,pad,mode='wrap')

    index[-pad:] += index[-pad-1]
    index[:pad]  -= index[-pad-1]

    for ii,idx in enumerate(index[pad:-pad]):

        # pool all values that falls within window
        pool = x[...,np.abs(index-idx)<=pad]

        std.append(np.full_like(pool[...,0],np.nanstd(pool)))

    return np.stack(std,axis=-1)

def running_regression_CV(x,y,z,index,window=30,lim=1,sub=np.nan):
    """
    Fits linear regression model: z = a * x + b * y, keeping dim -1 and -2.
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

    Returns 2 dimensional array (year,dayofyear)
    """
    slope_x,slope_y = [],[]

    pad   = window//2

    if len(x.shape)==2:
        x     = np.pad(x,pad,mode='wrap')[pad:-pad,:]
        y     = np.pad(y,pad,mode='wrap')[pad:-pad,:]
        z     = np.pad(z,pad,mode='wrap')[pad:-pad,:]

    if len(x.shape)==3:
        x     = np.pad(x,pad,mode='wrap')[pad:-pad,pad:-pad,:]
        y     = np.pad(y,pad,mode='wrap')[pad:-pad,pad:-pad,:]
        z     = np.pad(z,pad,mode='wrap')[pad:-pad,pad:-pad,:]


    index = np.pad(index,pad,mode='wrap')

    index[-pad:] += index[-pad-1]
    index[:pad]  -= index[-pad-1]

    for ii,idx in enumerate(index[pad:-pad]):

        # pool all values that falls within window
        xpool = x[...,np.abs(index-idx)<=pad]
        ypool = y[...,np.abs(index-idx)<=pad]
        zpool = z[...,np.abs(index-idx)<=pad]

        yslope_x,yslope_y = [],[]
        for yy in range(xpool.shape[-2]):

            # delete the relevant year from pool (for cross validation)
            filtered_xpool = np.delete(xpool,yy,axis=-2).flatten()
            filtered_ypool = np.delete(ypool,yy,axis=-2).flatten()
            filtered_zpool = np.delete(zpool,yy,axis=-2).flatten()

            idx_bool = np.logical_and(
                                np.logical_and(
                                                np.isfinite(filtered_xpool),
                                                np.isfinite(filtered_ypool)
                                              ),np.isfinite(filtered_zpool)
                            )

            if idx_bool.sum()<lim:
                yslope_x.append(sub)
                yslope_y.append(sub)

            else:
                X = np.stack(
                                [
                                    filtered_xpool[idx_bool],
                                    filtered_ypool[idx_bool]
                                ],axis=1
                            )

                Y = filtered_zpool[idx_bool]

                fit = sm.OLS(Y,X).fit()

                yslope_x.append(fit.params[0])
                yslope_y.append(fit.params[1])

        slope_x.append(np.array(yslope_x))
        slope_y.append(np.array(yslope_y))

    return np.stack(slope_x,axis=-1),np.stack(slope_y,axis=-1)

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


################################################################################
############################# Depricated #######################################
################################################################################
def persistence3(predictor,response,var,dim='time.dayofyear'):
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
