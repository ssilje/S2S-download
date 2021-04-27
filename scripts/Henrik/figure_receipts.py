# environment dependencies
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

# local dependencies
from S2S.local_configuration_H import config
from .create_domain_file import make_dir

def to_pandas(array,time,lon=[None],lat=[None]):
    """
    args:
        array: np.array, dim (step,time,...)
        time:  np.array, dim (step,time)

    returns

    """
    index   = ['step','time']
    df_base = {
                'step':np.stack([np.arange(1,array.shape[0]+1)]*array.shape[1],
                                                                        axis=1),
                'time':time,
                'values':array
            }

    if lon[0]:
        index.append('lon')
        lon = np.array(lon)
        df_base['step'] = np.stack([df_base['step']]*lon.shape[0],axis=-1)
        df_base['time'] = np.stack([df_base['time']]*lon.shape[0],axis=-1)
        df_base['lon']  = np.stack([np.stack([lon]*array.shape[1])]*array.shape[0])

    if lat[0]:
        index.append('lat')
        lat = np.array(lat)
        df_base['step'] = np.stack([df_base['step']]*lat.shape[0],axis=-1)
        df_base['time'] = np.stack([df_base['time']]*lat.shape[0],axis=-1)
        df_base['lon']  = np.stack([df_base['lon']]*lat.shape[0],axis=-1)
        df_base['lat']  = np.stack([np.stack([np.stack([lat]*lon.shape[0])]*\
                                                array.shape[1])]*array.shape[0])

    for key in df_base:
        df_base[key] = df_base[key].flatten()
    df_base['time'] = pd.to_datetime(df_base['time'])

    return pd.DataFrame.from_records(df_base,index=index)

def qq_plot(fc_anom,obs_anom,domainID):

    e_size   = fc_anom.shape[0]
    obs_anom = np.stack([obs_anom]*e_size,axis=0)

    fig,axes = plt.subplots(2,3,constrained_layout=True)

    for n,ax in enumerate(axes.flatten()):

        ax.scatter(y=fc_anom[:,n],x=obs_anom[:,n],alpha=0.6)

        ax.set_title('Lead time: '+str(n+1)+'W')
        ax.set_ylabel('Model')
        ax.set_xlabel('Observations')

    make_dir(config['SAVEFIG'])
    fig.suptitle(domainID)
    plt.savefig(config['SAVEFIG']+'test_qq.png')
    plt.close()

def crps_plot(score,time,domainID):

    fig,ax = plt.subplots(1,1,constrained_layout=True)
    ax.plot(time,score)
    plt.show()
    plt.close(0)
