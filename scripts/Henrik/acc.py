import numpy as np
import pandas as pd

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
