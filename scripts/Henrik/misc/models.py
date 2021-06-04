import numpy as np
import pandas as pd
import xarray as xr
import statsmodels.api as sm
import scipy.stats as stats
import time as time_lib

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
    res1 = []
    res2 = []
    res3 = []
    for n in range(y.shape[1]):
        X1 = np.nanmean(np.delete(x1,n,axis=1),axis=0).flatten()
        X2 = np.nanmean(np.delete(x2,n,axis=1),axis=0).flatten()
        X0 = np.ones_like(X1)
        X = np.stack(
                        [
                            X0,
                            X1,
                            X2
                        ],axis=1
                    )
        res = sm.OLS(np.nanmean(np.delete(y,n,axis=1),axis=0).flatten(), X).fit()
        res1.append(res.params[0])
        res2.append(res.params[1])
        res3.append(res.params[2])
    return np.array(res1),np.array(res2),np.array(res3)

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

def combo(observation,model,response,var,dim='time.dayofyear'):
    """
    Not very elegant and probably not particularly cheap
    """
    dim_name = dim.split('.')[0]
    print('\t performing models.combo()')
    observation = xr.concat([observation]*model.dims['member'],model.member)
    response    = xr.concat([response]*model.dims['member'],model.member)

    ds = xr.merge(
                    [
                        observation.rename({var:'observation'}),
                        model.rename({var:'model'}),
                        response.rename({var:'response'})
                ],join='inner',compat='override'
            )
    t = time_lib.time()
    N = len(list(ds.groupby(dim)))
    n = 0
    date_out,date_label_out = [],[]
    for date_label,date_group in list(ds.groupby(dim)):
        # n += 1
        # print('\t\t group ',n,' of ',N,
        #         ' total time: ',round(time_lib.time()-t,2))

        data = date_group
        tup = xr.apply_ufunc(
                regression2D, data.observation, data.model, data.response,
                input_core_dims  = [
                                    ['member',dim_name],
                                    ['member',dim_name],
                                    ['member',dim_name]
                                    ],
                output_core_dims = [[dim_name],[dim_name],[dim_name]],
                vectorize=True
                )

        date_out.append(xr.merge(
            [
                tup[0].rename('intercept'),
                tup[1].rename('slope_obs'),
                tup[2].rename('slope_model')
                ]
            )
        )
    return xr.concat(date_out,dim_name).sortby(dim_name)

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

################################################################################
############################# Depricated #######################################
################################################################################
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
