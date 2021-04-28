import pandas as pd
import numpy as np

def parse_month2season(months):

    season = {
    '1':'DJF',
    '2':'DJF',
    '3':'DJF',
    '4':'MAM',
    '5':'MAM',
    '6':'MAM',
    '7':'JJA',
    '8':'JJA',
    '9':'SON',
    '10':'SON',
    '11':'SON',
    '12':'DJF'
    }

    return [season[str(month)] for month in months]

def stack(array,dim,ax):
    return np.stack([array]*dim,axis=ax)

def fit_frame(array,position,shape):
    """
    args:
        array: np.array
        position: list containing position(s) of array in shape
        shape: list

    returns:
        array stacked as given by shape: np.array
    """
    position  = np.array(position)
    n_dims    = len(shape)

    for i in range(position[-1]+1,n_dims):
        array = stack(array,shape[i],-1)

    for i in range(position[0]-1,-1,-1):
        array = stack(array,shape[i],0)

    return array


def anom_to_pandas(fc_anom,obs_anom,time,lon,lat):
    """
    args:
        fc_anom:   np.array, dim (step,time,...)
        obs_anom:  np.array, dim (step,time)

    returns
        dataframe: pandas.DataFrame
    """

    shape = fc_anom.shape
    en,st,ti,lo,la = shape
    shape = list(shape)

    # stack observations to match ensemble
    obs_anom = stack(obs_anom,en,0)

    index   = ['ensemble_member','lead_time','time','lon','lat']

    df_base = dict.fromkeys(index)

    df_base['ensemble_member']  = fit_frame(np.arange(1,en+1),[0],shape)
    df_base['lead_time']        = fit_frame(np.arange(1,st+1),[1],shape)
    df_base['time']             = fit_frame(time[0],[2],shape)
    df_base['lon']              = fit_frame(lon,[3],shape)
    df_base['lat']              = fit_frame(lat,[4],shape)

    df_base['observations']     = obs_anom
    df_base['model']            = fc_anom

    for key in df_base:
        df_base[key] = df_base[key].flatten()
    df_base['time']  = pd.to_datetime(df_base['time'])
    df_base['month'] = df_base['time'].month_name()
    df_base['season'] = parse_month2season(df_base['time'].month)

    df = pd.DataFrame.from_records(df_base)

    return df

def score_to_pandas(array,time,lon=[None],lat=[None],name='values'):
    """
    args:
        array: np.array, dim (step,time,...)
        time:  np.array, dim (step,time)

    returns
        dataframe: pandas.DataFrame
    """
    index = ['lead_time','time']

    if lon[0]:
        index.append('lon')

    if lat[0]:
        index.append('lat')

    shape = list(array.shape)

    df_base = dict.fromkeys(index)

    df_base['lead_time']        = fit_frame(np.arange(1,shape[0]+1),[0],shape)
    df_base['time']             = fit_frame(time[0],[1],shape)

    if lon[0]:
        df_base['lon']          = fit_frame(lon,[2],shape)

    if lat[0]:
        df_base['lat']          = fit_frame(lat,[3],shape)

    df_base[name]               = array

    for key in df_base:
        df_base[key]  = df_base[key].flatten()
    df_base['time']   = pd.to_datetime(df_base['time'])
    df_base['month']  = df_base['time'].month_name()
    df_base['season'] = parse_month2season(df_base['time'].month)

    df = pd.DataFrame.from_records(df_base)

    return df
