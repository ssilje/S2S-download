import numpy as np
import pandas as pd
import xarray as xr
import statsmodels.api as sm

#
# Not cheap nor nice, TODO: re-do.
#

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

def persistence(predictor,response,var):
    """
    Not very elegant and probably not particularly cheap
    """
    ds = xr.merge(
                    [
                        predictor.rename({var:'predictor'}),
                        response.rename({var:'response'})
                ],join='inner',compat='override'
            )

    date_out,date_label_out = [],[]
    for date_label,date in list(ds.groupby('time.dayofyear')):
        loc_out,loc_label_out = [],[]
        for loc_label,loc in list(date.groupby('location')):
            step_out,step_label_out = [],[]
            for step_label,step in list(loc.groupby('step')):
                step_out.append(
                            a_and_b(
                                reg_m(
                                    step.response.values,
                                    [step.predictor.values]
                                )
                            )
                        )
                step_label_out.append(step_label)
            loc_out.append(xr.concat(step_out,pd.Index(step_label_out,name='step')))
            loc_label_out.append(loc_label)
        date_out.append(xr.concat(loc_out,pd.Index(loc_label_out,name='location')))
        date_label_out.append(date_label)
    return xr.concat(date_out,pd.Index(date_label_out,name='dayofyear'))

def combo(observation,model,response,var):
    """
    Not very elegant and probably not particularly cheap
    """
    observation = xr.concat([observation]*model.dims['member'],model.member)
    response    = xr.concat([response]*model.dims['member'],model.member)

    ds = xr.merge(
                    [
                        observation.rename({var:'observation'}),
                        model.rename({var:'model'}),
                        response.rename({var:'response'})
                ],join='inner',compat='override'
            )

    date_out,date_label_out = [],[]
    for date_label,date in list(ds.groupby('time.dayofyear')):
        loc_out,loc_label_out = [],[]
        for loc_label,loc in list(date.groupby('location')):
            step_out,step_label_out = [],[]
            for step_label,step in list(loc.groupby('step')):
                step_out.append(
                            a_b_and_c(
                                reg_m(
                                    step.response.values.flatten(),
                                    [
                                        step.observation.values.flatten(),
                                        step.model.values.flatten()
                                    ]
                                )
                            )
                        )
                step_label_out.append(step_label)
            loc_out.append(xr.concat(step_out,pd.Index(step_label_out,name='step')))
            loc_label_out.append(loc_label)
        date_out.append(xr.concat(loc_out,pd.Index(loc_label_out,name='location')))
        date_label_out.append(date_label)
    return xr.concat(date_out,pd.Index(date_label_out,name='dayofyear'))
