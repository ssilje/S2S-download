import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
import pandas as pd

import scripts.Henrik.latex as latex
from S2S.local_configuration import config
import scripts.Henrik.xarray_helpers as xh
from scripts.Henrik.organize_barentzwatch import name_from_loc

def make_dir(path):
    """
    Creates directory if it does not exist.

    args:
        path: str
    """
    if not os.path.exists(path):
        os.makedirs(path)

def save_fig(fig,filename):
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
    M = np.array(m)
    while n > M[0]*M[1]:
        if M[0]==M[1]:
            M[1] = M[1]+1
        else:
            M[0] = M[0]+1
    return tuple(M)

def fg(da,dim):
    return n2m(len(list(da.groupby(dim))))

def qq_plot(dax,day,dim='validation_time.month',
                        x_axlabel='obs',
                        y_axlabel='model',
                        filename='qq_plot',
                        title=''):

    dax,day  = xr.broadcast(dax,day)

    if dim=='validation_time.month':
        dax = xh.assign_validation_time(dax)
        day = xh.assign_validation_time(day)

    filename = filename+'_'+name_from_loc(str(dax.location.values))
    title    = title + ' ' + \
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

        ax = axes[n]

        ax.set_title(month(xlabel))

        ax.set_xlabel(x_axlabel)
        ax.set_ylabel(y_axlabel)

        ax.plot(x,y,'o',alpha=0.4,ms=1)
        ax.plot([0, 1], [0, 1],'k',transform=ax.transAxes,alpha=0.7,linewidth=0.6)
        ax.axis('equal')

    fig.suptitle(title)
    save_fig(fig,filename)

def skill_plot(mod,clim,dim='validation_time.month',
                                filename='crpss',
                                title=''):

    if dim=='validation_time.month':
        mod = xh.assign_validation_time(mod)
        clim = xh.assign_validation_time(clim)

    filename = filename+'_'+name_from_loc(str(clim.location.values))
    title    = title+' '+name_from_loc(str(clim.location.values))

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
        z = ss.values[0]

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

    fig.suptitle(title)
    save_fig(fig,filename)
