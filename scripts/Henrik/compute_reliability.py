import numpy as np
import xarray as xr
import pandas as pd
from xskillscore import reliability

def exceedance(dataarray,threshold,direction=1):

    if direction < 0:
        dataarray = dataarray.where(dataarray < threshold,np.nan)
    if direction > 0:
        dataarray = dataarray.where(dataarray > threshold,np.nan)

    return dataarray.where(np.isnan(dataarray),np.float32(1.))\
                    .where(np.logical_not(np.isnan(dataarray)),np.float32(0.))

def bool_exceedance(dataarray,threshold,direction=1):

    if direction < 0:
        return dataarray < threshold
    if direction > 0:
        return dataarray > threshold

def probability_of_exceedance(dataarray,threshold,direction=1):

    dataarray = exceedance(dataarray,threshold,direction)

    return dataarray.mean('ensemble_member')

def reliability_of_forecast(
                            forecast,
                            observations,
                            thresholds,
                            direction=1,
                            dim='time',
                            probability_bin_edges=np.arange(0,1.1,0.1)
                            ):

    index = pd.Index(thresholds,name='threshold')
    out = []
    for threshold in thresholds:
        obs    = bool_exceedance(observations,threshold,threshold)
        fc_prb = bool_exceedance(forecast,threshold,threshold)\
                                                .mean('ensemble_member')

        out.append(
                reliability(
                            obs,
                            fc_prb,
                            dim,
                            probability_bin_edges
                        )
                    )

    return xr.concat(out,index)

def grouped_reliability_of_forecast(
                                    forecast,
                                    observations,
                                    thresholds,
                                    direction=1,
                                    dim='time',
                                    group='time.month',
                                    probability_bin_edges=np.arange(0,1.2,0.2)
                                ):
    out_step     = []
    index_step   = []
    grouped_obs_step = list(observations.groupby('step'))
    grouped_fc_step  = list(forecast.groupby('step'))
    for s in range(len(grouped_obs_step)):

        index_step.append(grouped_fc_step[s][0])

        out_group   = []
        index_group = []
        grouped_obs = list(grouped_obs_step[s][1].groupby(group))
        grouped_fc  = list(grouped_fc_step[s][1].groupby(group))
        for n in range(len(grouped_obs)):

            # change such that season/threshold is arbitrary
            label = grouped_fc[n][0]
            index_group.append(label)
            out_group.append(
                    reliability_of_forecast(
                                    grouped_fc[n][1],
                                    grouped_obs[n][1],
                                    np.array(thresholds.sel(season=label)),
                                    direction,
                                    dim=dim,
                                    probability_bin_edges=probability_bin_edges
                                    )
                                )
        out_step.append(out_group)
    return out_step,index_group,index_step

def compute(forecast,observations,threshold):
    """
    args:
        forecast:       np.array, dim (ensemble,step,time,...)
        observation:    np.array, dim (step,time,...)
        threshold:      np.array, dim (n,step,time,...)

    returns:
        reliability:    np.array, dim (n,...,2)

    - Dimensions (step,time,...) must be identical.
    """
    p_bins            = np.arange(0,1.1,0.1)
    ensemble_size     = forecast.shape[0]
    f_times           = forecast.shape[2]

    if threshold.shape == observations.shape:
        threshold = np.array([threshold])

    if not threshold[0].shape == observations.shape:
        print('dimension mismatch')
        exit()

    out = []
    for t in threshold:

        if t.flatten()[0] < 0:
            hits_fc_ens  = forecast     <= t
            hits_obs     = observations <= t
        else:
            hits_fc_ens  = forecast     >= t
            hits_obs     = observations >= t

        prob_fc  = hits_fc_ens.sum(axis=0)/ensemble_size

        obs_freq          = []
        for p in range(len(p_bins)-1):

            fc_hit_idx = np.ones_like(observations,dtype=int)

            fc_hit_idx[prob_fc >= p_bins[p+1]] = 0
            fc_hit_idx[prob_fc < p_bins[p]] = 0
            print(fc_hit_idx.shape)
            print(hits_obs.shape)
            obs_freq.append(
                            hits_obs[fc_hit_idx].sum(axis=1)/\
                            fc_hit_idx.sum(axis=1)
                            )
        out.append(np.array(obs_freq))
    return np.array(out),p_bins
