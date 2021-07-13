import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
import pandas as pd
import xskillscore as xs
import json

from S2S.local_configuration import config

import .latex as latex
import scripts.Henrik.xarray_helpers as xh

def name_from_loc(loc):
    """
    Returns name of Barentswatch location from loc number
    """
    with open(config['SITES'], 'r') as file:
        data = json.load(file)
        for line in data:
            if line["localityNo"]==int(loc):
                return line['name']

def make_dir(path):
    """
    Creates directory if it does not exist.

    args:
        path: str
    """
    if not os.path.exists(path):
        os.makedirs(path)

def save_fig(fig,filename):
    """
    Save figure as pdf and png
    """
    make_dir(config['SAVEFIG'])
    # file,converge = latex.save_figure(fig,config['SAVEFIG']+filename)
    plt.savefig(config['SAVEFIG']+\
                    filename+'.png',dpi='figure',bbox_inches='tight')
    plt.close()
    print('Figure stored at: '+config['SAVEFIG']+filename+'.pdf')

def month(ii):
    """
    ii : integer

    returns

    month string_tools
    """
    ii = int(ii)
    return ['JAN','FEB','MAR',\
            'APR','MAI','JUN',\
            'JUL','AUG','SEP',\
            'OCT','NOV','DES'][ii-1]


def n2m(n,m=(1,1)):
    """
    Returns minimum required arangement of subplots
    for n differents axis
    """
    M = np.array(m)
    while n > M[0]*M[1]:
        if M[0]==M[1]:
            M[1] = M[1]+1
        else:
            M[0] = M[0]+1
    return tuple(M)

def fg(da,dim):
    """
    Wrapper for n2m()
    """
    return n2m(len(list(da.groupby(dim))))

def qq_plot(dax_in,day_in,dim='validation_time.month',
                        x_axlabel='obs',
                        y_axlabel='model',
                        filename='',
                        title=''):

    for loc in dax_in.location:

        dax = dax_in.sel(location=loc)
        day = day_in.sel(location=loc)

        dax,day  = xr.broadcast(dax,day)

        if dim=='validation_time.month':
            dax = xh.assign_validation_time(dax)
            day = xh.assign_validation_time(day)

        fname    = 'qqplot_'+filename+'_'+name_from_loc(str(dax.location.values))
        suptitle = title + ' ' + \
                        name_from_loc(str(dax.location.values)) #+\
                        # ' leadtime: ' + '-'.join([
                        #         str(pd.Timedelta(dax.step.min().to_pandas()).days),
                        #         str(pd.Timedelta(dax.step.max().to_pandas()).days)
                        #     ]) + ' days'
        ###########################
        #### Initialize figure ####
        ###########################
        latex.set_style(style='white')
        subplots = fg(dax,dim)
        fig,axes = plt.subplots(subplots[0],subplots[1],
                        figsize=latex.set_size(width='thesis',
                            subplots=(subplots[0],subplots[1]))
                        )
        axes = axes.flatten()
        ###########################

        x_group = list(dax.groupby(dim))
        y_group = list(day.groupby(dim))

        for n,(xlabel,xdata) in enumerate(x_group):

            ylabel,ydata = y_group[n]
            # this approach is not bulletproof
            # check manually that right groups go together
            print('\t graphics.qq_plot: matched groups ',xlabel,ylabel)

            x,y = xr.align(xdata,ydata)
            x,y = x.values.flatten(),y.values.flatten()
            all = np.array([x,y]).flatten()

            ax = axes[n]

            ax.set_title(month(xlabel))

            ax.set_xlabel(x_axlabel)
            ax.set_ylabel(y_axlabel)

            limits = [np.nanmin(all)-1,np.nanmax(all)+1]
            ax.set_aspect('equal')
            ax.set_xlim(limits)
            ax.set_ylim(limits)
            ax.plot([0, 1], [0, 1],'k',transform=ax.transAxes,alpha=0.7,linewidth=0.6)

            ax.plot(x,y,'o',alpha=0.4,ms=1)

        fig.suptitle(suptitle)
        save_fig(fig,fname)

def skill_plot(in_mod,in_clim,dim='validation_time.month',
                                filename='',
                                title=''):
    for loc in in_mod.location:

        mod  = in_mod.sel(location=loc)
        clim = in_clim.sel(location=loc)

        if dim=='validation_time.month':
            mod = xh.assign_validation_time(mod)
            clim = xh.assign_validation_time(clim)

        fname    = 'crpss'+\
                        '_'+filename+'_'+\
                                name_from_loc(str(clim.location.values))
        suptitle = title+' '+name_from_loc(str(clim.location.values))

        dim_name = dim.split('.')[0]
        ###########################
        #### Initialize figure ####
        ###########################
        latex.set_style(style='white')
        subplots = fg(mod,dim)
        fig,axes = plt.subplots(subplots[0],subplots[1],
                        figsize=latex.set_size(width='thesis',
                            subplots=(subplots[0],subplots[1]))
                        )
        axes = axes.flatten()
        ###########################

        x_group = list(mod.groupby(dim))
        y_group = list(clim.groupby(dim))

        for n,(xlabel,xdata) in enumerate(x_group):

            ylabel,ydata = y_group[n]
            # this approach is not bulletproof
            # check manually that right groups go together
            print('\t graphics.skill_plot: matched groups ',xlabel,ylabel)

            xdata = xdata.unstack().sortby(['time','step'])
            ydata = ydata.unstack().sortby(['time','step'])

            xdata,ydata = xr.align(xdata,ydata)

            ss = 1 - xdata/ydata
            ss = ss.mean('time',skipna=True)

            x = np.array([td.days for td in ss.step.to_pandas()])
            z = ss.values

            ax = axes[n]

            ax.set_title(month(xlabel))
            ax.set_xlabel('lead time [D]')
            ax.set_ylabel('CRPSS')

            ax.plot(x,z,'o-',alpha=0.4,ms=1)
            # ax.plot(x,zx,'o-',alpha=0.4,ms=1,label='x')
            # ax.plot(y,zy,'o-',alpha=0.4,ms=1,label='y')
            # ax.legend()
            ax.set_ylim((-1,1))
            ax.plot([0, 1], [0.5, 0.5],'k',transform=ax.transAxes,alpha=0.7,linewidth=0.6)

        fig.suptitle(suptitle)
        save_fig(fig,fname)

def timeseries(
            observations,
            cast,
            lead_time=[7,14,21],
            olab='obs',
            clabs=['forecast'],
            filename='',
            title=''
            ):



    for loc in observations.location:

        latex.set_style(style='white')
        fig,axes = plt.subplots(len(lead_time),1,
                        figsize=latex.set_size(width='thesis',
                            subplots=(1,1))
                        )
        clabs = clabs*len(cast)

        fname = 'timeseries_'+filename+'_'+str(name_from_loc(loc.values))

        for n,lt in enumerate(lead_time):

            fname += '_'+str(lt)

            if hasattr(axes,'__iter__'):
                ax = axes[n]
            else:
                ax = axes


            ax.set_xlabel('Forecasted time')
            ax.set_ylabel('')

            o   = observations\
                    .sel(step=pd.Timedelta(lt,'D'),location=loc)\
                        .sortby('time')

            ax.plot(
                    o.validation_time.squeeze(),
                    o.squeeze(),
                    'o',
                    color='k',
                    ms=0.2,
                    linewidth=0.5,
                    label=olab
                )

            subtitle = '    MAE:'
            for cn,c in enumerate(cast):

                c = c\
                        .sel(step=pd.Timedelta(lt,'D'),location=loc)\
                            .sortby('time')
                c = c.sel(time=slice(o.time.min(),o.time.max()))

                ax.plot(
                        c.validation_time.squeeze(),
                        c.mean('member').squeeze(),
                        alpha=0.7,
                        label=clabs[cn],
                        linewidth=0.5,
                        zorder=30
                    )
                std = c.std('member').squeeze()
                ax.fill_between(
                        c.validation_time.squeeze(),
                        c.mean('member').squeeze()+std,
                        c.mean('member').squeeze()-std,
                        alpha=0.3,
                        zorder=30
                    )

                score = xs.mae(o,c.mean('member'),dim=[])\
                            .mean('time',skipna=True)

                subtitle += ' '+clabs[cn]+' '+str(round(float(score.values),2))

            ax.set_title('Lead time '+str(lt)+' days '+subtitle)
            ax.legend()

        fig.suptitle(str(name_from_loc(loc.values))+' '+title)
        save_fig(fig,fname)
