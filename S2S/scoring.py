import numpy as np
import pandas as pd
import xarray as xr

def ACC(fc_anom,obs_anom,weights):
    """
    Anomaly correlation after Wilks (2011, Chapter 8)
    D.S. Wilks, Chapter 8 - Forecast Verification,
    International Geophysics, Academic Press, Volume 100, 2011,
    Pages 301-394, https://doi.org/10.1016/B978-0-12-385022-5.00008-7.
    args:
        fc_anom:  np.array, dims (ensemble_member,step,validation_time,lon,lat)
                              or (step,validation_time,lon,lat)
        obs_anom: np.array, dims (step,validation_time,lon,lat)
    returns:
        acc:      np.array, dims (step,validation_time)
    """
    # mean over ensemble members
    if fc_anom.shape[0] != obs_anom.shape[0]:
        fc_anom = fc_anom.mean(axis=0)

    if weights.shape[0] != obs_anom.shape[0]:
        weights = weights.mean(axis=0)

    if fc_anom.shape != obs_anom.shape:
        print('Input dim mismatch')
        exit()

    not_nan = np.isfinite(fc_anom)

    sum_weights = (weights*not_nan).sum(axis=(-1,-2))

    fc_mean     = np.nansum(fc_anom*not_nan*weights,axis=(-1,-2))/sum_weights
    obs_mean    = np.nansum(obs_anom*not_nan*weights,axis=(-1,-2))/sum_weights

    fc_anom  = np.transpose(fc_anom)
    obs_anom = np.transpose(obs_anom)

    fc_mean  = np.transpose(fc_mean)
    obs_mean = np.transpose(obs_mean)

    var_fc  = (fc_anom-fc_mean)
    var_obs = (obs_anom-obs_mean)

    covar   = var_fc*var_obs
    not_nan = np.isfinite(covar)

    return\
        np.transpose(\
            np.nansum(covar*not_nan,axis=(0,1))/\
                np.sqrt(
                    np.nansum((var_fc*not_nan)**2,axis=(0,1))*\
                    np.nansum((var_obs*not_nan)**2,axis=(0,1)\
                )
            )
        )




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

    obs = obs.broadcast_like(fc.mean('member'))

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

def fair_factor(ensemble_size):
    return ensemble_size/(ensemble_size+1)
