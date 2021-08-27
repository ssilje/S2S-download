import numpy as np
import pandas as pd
import xarray as xr
import xskillscore as xs

import S2S.xarray_helpers as xh

def centered_acc(x,y):
    """
    Anomaly correlation after Wilks (2011, Chapter 8) - Eq. 8.

    D.S. Wilks, Chapter 8 - Forecast Verification,
    International Geophysics, Academic Press, Volume 100, 2011,
    Pages 301-394, https://doi.org/10.1016/B978-0-12-385022-5.00008-7.
    """
    idx_bool = ~np.logical_or(np.isnan(x),np.isnan(y))

    x = x[idx_bool]
    y = y[idx_bool]

    M = len(x)

    x_anom = (x - x.mean())
    y_anom = (y - y.mean())

    covar = x_anom * y_anom / M

    var_x = x_anom**2 / M
    var_y = y_anom**2 / M

    return covar.sum() / np.sqrt( ( var_x.sum() + var_y.sum() ) )

def uncentered_acc(x,y):
    """
    Anomaly correlation after Wilks (2011, Chapter 8) - Eq. 8.64

    D.S. Wilks, Chapter 8 - Forecast Verification,
    International Geophysics, Academic Press, Volume 100, 2011,
    Pages 301-394, https://doi.org/10.1016/B978-0-12-385022-5.00008-7.
    """
    idx_bool = ~np.logical_or(np.isnan(x),np.isnan(y))

    x = x[idx_bool]
    y = y[idx_bool]

    M = len(x)

    x_anom = x
    y_anom = y

    covar = x_anom * y_anom / M

    var_x = x_anom**2 / M
    var_y = y_anom**2 / M

    return covar.sum() / np.sqrt( ( var_x.sum() + var_y.sum() ) )

def ACc(forecast,observations,weights=None,centered=True):
    """
    Anomaly correlation after Wilks (2011, Chapter 8)

    D.S. Wilks, Chapter 8 - Forecast Verification,
    International Geophysics, Academic Press, Volume 100, 2011,
    Pages 301-394, https://doi.org/10.1016/B978-0-12-385022-5.00008-7.
    """
    if weights is None:
        weights = xr.full_like(observations,1.)

    forecast     = forecast * weights
    observations = observations * weights

    try:
        forecast = forecast.mean('member')
    except AttributeError:
        pass


    ds = xr.merge(
                    [
                        forecast.rename('fc'),
                        observations.rename('obs')
                    ],join='inner',compat='override'
                )

    if centered:
        ufunc = centered_acc
    else:
        ufunc = uncentered_acc

    r = xr.apply_ufunc(
            ufunc,ds.fc,ds.obs,
            input_core_dims = [['lat','lon'],['lat','lon']],
            output_core_dims = [[]],
            vectorize=True,dask='parallelized'
        )

    return r



def CRPS_ensemble(obs,fc,fair=True,axis=0):
    """
    @author: Ole Wulff
    @date: 2020-07-08

    implementation of fair (adjusted) CRPS based on equation (6) from Leutbecher (2018, QJRMS, https://doi.org/10.1002/qj.3387)
    version with fair=False tested against properscoring implementation crps_ensemble (see https://pypi.org/project/properscoring/)

    INPUT:
        obs: observations as n-dimensional array
        fc: forecast ensemble as (n+1)-dimensional where the extra dimension (axis) carries the ensemble members
        fair: if True returns the fair version of the CRPS accounting for the limited ensemble size (see Leutbecher, 2018)
              if False returns the normal CRPS
        axis: axis of fc array that contains the ensemble members, defaults to 0
    OUTPUT:
        CRPS: n-dimensional array
    TODO:
        implement weights for ensemble member weighting
    """
    odims = obs.shape
    M = fc.shape[axis]
    if axis != 0:
        fc = np.swapaxes(fc,axis,0)

    # flatten all dimensions except for the ensemble member dimension:
    fc_flat = fc.reshape([M,-1])
    obs_flat = obs.reshape([-1])

    dsum = np.array([abs(fc_flat[jj] - fc_flat[kk]) for kk in range(M) for jj in range(M)]).sum(axis=axis)
    if fair:
        CRPS = 1/M * (abs(fc_flat - obs_flat)).sum(axis=axis) - 1/(2*M*(M-1)) * dsum
    else:
        CRPS = 1/M * (abs(fc_flat - obs_flat)).sum(axis=axis) - 1/(2*M**2) * dsum

    # is this necessary or even a good idea at all?
#     del dsum, fc_flat, obs_flat

    return CRPS.reshape([*odims])

def crps_ensemble(obs,fc,fair=True):
    """
    A xarray wrapper for CRPS_ensemble()
    """

    # obs = obs.broadcast_like(fc.mean('member'))
    obs,fc = xr.align(obs,fc,join='outer')
    print(obs)

    return xr.apply_ufunc(
            CRPS_ensemble, obs, fc,
            input_core_dims  = [['time'],['member','time']],
            output_core_dims = [['time']],
            vectorize=True
        )

def fair_brier_score(observations,forecasts):

    M = forecasts['member'].size
    e = (forecasts == 1).sum('member')
    o = observations
    return (e / M - o) ** 2 - e * (M - e) / (M ** 2 * (M - 1))

class SSCORE:
    """
    Calculate monthly and seasonal means of scores with bootstrapped CIs.
    """

    def __init__(self,**kwargs):

        self.forecast     = self.try_key(kwargs,'forecast')
        self.observations = self.try_key(kwargs,'observations')

        self.data  = xr.merge(
                            [
                                self.forecast.rename('fc'),
                                self.observations.rename('obs')
                            ],
                        join='inner',
                        compat='equals'
                        )

    @staticmethod
    def skill_score(x,y):
        return 1 - np.nanmean(x,axis=-1)/np.nanmean(y,axis=-1)

    @staticmethod
    def try_key(dictionary,key):
        try:
            return dictionary[key]
        except KeyError:
            return None

    def bootstrap(self,N=1000,ci=.95,min_period=2):

        # split into weeks
        ds = xh.unstack_time(self.data)

        low_q,est,high_q,ny = xr.apply_ufunc(
                self.pull, ds.fc, ds.obs, N, ci, min_period,
                input_core_dims  = [
                                    ['dayofyear','year'],
                                    ['dayofyear','year'],
                                    [],
                                    [],
                                    []
                                ],
                output_core_dims = [[],[],[],[]],
                vectorize=True
            )
        out = xr.merge(
                    [
                        low_q.rename('low_q'),
                        est.rename('est'),
                        high_q.rename('high_q')
                    ], join='inner', compat='equals'
                )
        out = out.assign_coords(number_of_years=ny)
        return out

    def pull(self,fc,obs,N,ci=.95,min_period=2):

        # use only years of at least min_period observations
        fc_idx  = (np.isfinite(fc)).sum(axis=0)>min_period
        obs_idx = (np.isfinite(obs)).sum(axis=0)>min_period

        # use only coincinding years
        idx_bool = np.logical_and(
                            fc_idx,
                            obs_idx
                            )

        ny = idx_bool.sum()

        if idx_bool.sum()>1:

            fc = np.nanmean(fc[...,idx_bool],axis=0)
            obs = np.nanmean(obs[...,idx_bool],axis=0)

            y  = fc.shape[-1]

            # generate random integers as indices for fc-obs pairs
            _idx_ = np.random.randint(low=0,high=y,size=(N,y))

            # pick y random fc-obs pairs N times
            _fc_   = fc[_idx_]
            _obs_  = obs[_idx_]

            # calculate score N times
            score = np.sort(self.skill_score(_fc_,_obs_))

            # quantiles
            alpha  = (1-ci)/2

            high_q = score[ int( (N-1) * (1-alpha) ) ]
            low_q  = score[ int( (N-1) * alpha ) ]

        else:

            high_q = np.nan
            low_q  = np.nan

        if idx_bool.sum()>0:
            # actual score
            est_score = self.skill_score(fc,obs)
        else:
            est_score = np.nan

        return low_q,est_score,high_q,ny
