from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
import pandas as pd
import xskillscore as xs
import json
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl

from S2S.local_configuration import config
import S2S.xarray_helpers as xh
from .graphics import *
import S2S.location_cluster as lc
from S2S.scoring import SSCORE

from . import latex

def ss(fc,clim,time='before',est='mean'):

    if time=='before':

        if est=='mean':
            fc   = fc.mean('time',skipna=True)
            clim = clim.mean('time',skipna=True)
        elif est=='median':
            fc   = fc.median('time',skipna=True)
            clim = clim.median('time',skipna=True)

    SS = 1 - fc/clim

    if time=='after':

        if est=='mean':
            SS = SS.mean('time',skipna=True)
        elif est=='median':
            SS = SS.median('time',skipna=True)

    return SS

def skill_agg(
            in_obs,
            in_mod,
            clim_mean,
            dim='validation_time.month',
            filename='',
            title='',
            ylab='',
            mlabs=[''],
            mcols=['blue','orange','green','red'],
            aggregate_around=None,
        ):

    mcols = mcols*len(in_mod)
    mlabs = mlabs*len(in_mod)

    # if location dimension does not exist, assign it
    if np.isin(np.array(in_obs.dims),'location').sum() == 0:

        out = []
        for array in in_mod:

            out.append(array.expand_dims('location'))

        clim_mean = clim_mean.expand_dims('location')
        in_obs    = in_obs.expand_dims('location')

        in_mod = out

    if aggregate_around is not None:

        locations = [aggregate_around]

    else:

        if in_obs.location.values.ndim==0:
            locations = [in_obs.location]
        else:
            locations = in_obs.location

    for loc in locations:

        if aggregate_around is None:
            clim = xh.assign_validation_time(in_obs).sel(location=loc)
        else:
            clim = xh.assign_validation_time(in_obs)

        fname    = 'mae/'+filename+'_'+dim.split('.')[1]+'_'+\
                                name_from_loc(str(loc.values))
        suptitle = title+' '+name_from_loc(str(loc.values))

        dim_name = dim.split('.')[0]
        ###########################
        #### Initialize figure ####
        ###########################
        latex.set_style(style='white')
        subplots = fg(clim,dim)
        fig,axes = plt.subplots(subplots[0],subplots[1],
                        figsize=latex.set_size(width='thesis',
                            subplots=(subplots[0],subplots[1]))
                        )
        axes = axes.flatten()
        ###########################

        for model,mlab,mcol in zip(in_mod,mlabs,mcols):

            if aggregate_around is None:
                mod = xh.assign_validation_time(model.sel(location=loc))
                cm  = xh.assign_validation_time(clim_mean.sel(location=loc))

            else:
                mod = xh.assign_validation_time(model)
                cm  = xh.assign_validation_time(clim_mean)


            x_group  = list(mod.groupby(dim))
            y_group  = list(clim.groupby(dim))
            cm_group = list(cm.groupby(dim))

            for n,(xlabel,xdata) in enumerate(x_group):

                ax = axes[n]

                ylabel,ydata   = y_group[n]
                cmlabel,cmdata = cm_group[n]
                # this approach is not bulletproof
                # check manually that right groups go together
                print('\tgraphics.skill_plot: matched groups ',
                                            xlabel,ylabel,cmlabel)

                xdata  = xdata.unstack().sortby(['time','step'])
                ydata  = ydata.unstack().sortby(['time','step'])
                cmdata = cmdata.unstack().sortby(['time','step'])

                xdata,ydata,cmdata = xr.align(xdata,ydata,cmdata)

                lead_time = np.array([td.days for td in ydata.step.to_pandas()])

                try:
                    ens_center_estimate = xdata.mean('member',skipna=True)
                except ValueError:
                    ens_center_estimate = xdata

                try:
                    ens_center_estimate_cm = cmdata.mean('member',skipna=True)
                except ValueError:
                    ens_center_estimate_cm = cmdata

                score_fc   = xs.mae(
                                    ens_center_estimate,
                                    ydata,
                                    dim=[]
                                )
                score_clim = xs.mae(ens_center_estimate_cm,ydata,dim=[])

                if aggregate_around is not None:

                    score_fc = score_fc.mean('location',skipna=True)
                    score_clim = score_clim.mean('location',skipna=True)

                if dim.split('.')[1]=='season':
                    min_period=6
                else:
                    min_period=2

                SS = SSCORE(
                            observations = score_clim,
                            forecast     = score_fc
                        ).bootstrap( N = 10000, min_period = min_period )

                ax.plot(
                        lead_time,
                        SS.low_q,
                        '--',
                        color=mcol,
                        linewidth=0.7,
                        alpha=0.7,
                        label='95% CI'
                        )
                ax.plot(
                        lead_time,
                        SS.high_q,
                        '--',
                        color=mcol,
                        linewidth=0.7,
                        alpha=0.7,
                        label='95% CI'
                        )
                ax.plot(
                        lead_time,
                        SS.est,
                        '-',
                        color=mcol,
                        linewidth=0.9,
                        alpha=0.9,
                        label='MAE SS est.' + mlab
                        )

                ax.fill_between(
                        lead_time,
                        SS.high_q,
                        SS.low_q,
                        alpha=0.1,
                        zorder=30
                        )

                for lt in lead_time:
                    ax.text(
                            lt,1.1,#0.83,
                            str(
                                SS.number_of_years\
                                    .sel(step=pd.Timedelta(lt,'D')).values
                            ),
                            horizontalalignment='center',
                            verticalalignment='center',
                            color='red',label='number of years for est.'
                        )
                plt.plot([], [], ' ',
                            label='Red: Number of years per estimate'
                        )

                ax.set_ylim((-1,1))
                ax.set_xlim((lead_time[0]-1,lead_time[-1]+1))

                ax.set_title(month(xlabel))
                ax.set_ylabel(ylab)
                ax.set_xticks(lead_time)

                ax.plot(
                        [0, 1],[0.5, 0.5],'k',
                        transform=ax.transAxes,alpha=0.7,
                        linewidth=0.6
                        )

                row = n//subplots[1]
                col = n%subplots[1]

                if col==0:
                    ax.set_ylabel('MAE SS')
                if row==subplots[0]-1:
                    ax.set_xlabel('lead time [D]')


        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        legend = fig.legend(
                    by_label.values(),
                    by_label.keys(),
                    loc='upper center',
                    bbox_to_anchor=(0.5, -0.05),
                    fancybox=True,
                    shadow=True,
                    ncol=4
                )
        legend.get_texts()[0].set_color('r')

        fig.suptitle(suptitle)
        save_fig(fig,fname)

def skill_agg_cluster(
            in_obs,
            in_mod,
            clim_mean,
            clim_std,
            loc,
            dim='validation_time.month',
            filename='',
            title='',
            ylab=''
        ):

    fname    = 'mae_cluster/'+filename+'_'+loc+'_'+dim.split('.')[1]
    suptitle = title+' '+loc

    dim_name = dim.split('.')[0]
    ###########################
    #### Initialize figure ####
    ###########################
    latex.set_style(style='white')
    subplots = fg(in_obs,dim)
    fig,axes = plt.subplots(subplots[0],subplots[1],
                    figsize=latex.set_size(width='thesis',
                        subplots=(subplots[0],subplots[1]))
                    )
    axes = axes.flatten()
    ###########################

    for model in in_mod:

        mod  = model
        cm   = clim_mean
        cs   = clim_std
        clim = in_obs

        x_group  = list(mod.groupby(dim))
        y_group  = list(clim.groupby(dim))
        cm_group = list(cm.groupby(dim))
        cs_group = list(cs.groupby(dim))

        for n,(xlabel,xdata) in enumerate(x_group):

            ax = axes[n]

            ylabel,ydata   = y_group[n]
            cmlabel,cmdata = cm_group[n]
            cslabel,csdata = cs_group[n]
            # this approach is not bulletproof
            # check manually that right groups go together
            print('\tgraphics.skill_plot: matched groups ',
                                        xlabel,ylabel,cmlabel,cslabel)

            xdata  = xdata.unstack().sortby(['time','step'])
            ydata  = ydata.unstack().sortby(['time','step'])
            cmdata = cmdata.unstack().sortby(['time','step'])
            csdata = csdata.unstack().sortby(['time','step'])

            xdata,ydata,cmdata,csdata = xr.align(xdata,ydata,cmdata,csdata)

            lead_time = np.array([td.days for td in ydata.step.to_pandas()])

            # score_cf     = xs.mae(xdata.sel(member=0),ydata,dim=[])
            score_mean = xs.mae(xdata.mean('member',skipna=True),ydata,dim=[])
            # score_median = xs.mae(xdata.median('member',skipna=True),ydata,dim=[])

            score_clim = xs.mae(cmdata,ydata,dim=[])

            ax.plot(
                    lead_time,
                    ss(
                        score_mean.mean('location',skipna=True),
                        score_clim.mean('location',skipna=True),
                        time='before',
                        est='mean'
                    ),
                    '--',color='k',linewidth=0.7,
                    alpha=0.5,
                    label='before mean'
                    )
            ax.plot(
                    lead_time,
                    ss(
                        score_mean,
                        score_clim,
                        time='after',
                        est='mean'
                    ).mean('location',skipna=True),
                    '-',color='k',linewidth=0.7,
                    alpha=0.5,
                    label='after mean'
                    )
            ax.plot(
                    lead_time,
                    ss(
                        score_mean.mean('location',skipna=True),
                        score_clim.mean('location',skipna=True),
                        time='before',
                        est='median'
                    ),
                    '--',color='r',linewidth=0.7,
                    alpha=0.5,
                    label='before median'
                    )
            ax.plot(
                    lead_time,
                    ss(
                        score_mean,
                        score_clim,
                        time='after',
                        est='median'
                    ).mean('location',skipna=True),
                    '-',color='r',linewidth=0.7,
                    alpha=0.5,
                    label='after median'
                    )

            ax.set_ylim((-1,1))

            ax.set_title(month(xlabel))
            ax.set_ylabel(ylab)
            ax.set_xticks(lead_time)

            ax.plot(
                    [0, 1],[0.5, 0.5],'k',
                    transform=ax.transAxes,alpha=0.7,
                    linewidth=0.6
                    )

            row = n//subplots[1]
            col = n%subplots[1]

            if col==0:
                ax.set_ylabel('MAE SS')
            if row==subplots[0]-1:
                ax.set_xlabel('lead time [D]')


    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    fig.legend(
                by_label.values(),
                by_label.keys(),
                loc='upper center',
                bbox_to_anchor=(0.5, -0.05),
                fancybox=True,
                shadow=True,
                ncol=4
            )

    fig.suptitle(suptitle)
    save_fig(fig,fname)

def map(
              observations,
              model,
              clim_mean,
              dim='validation_time.month',
              title='',
              lead_time=[9,16,23,30,37],
              filename=''
             ):

    for lt in lead_time:

        fname = 'mae_map/'+filename+'_'+dim.split('.')[1]+'_'+str(lt)

        mod = model.sel(step=pd.Timedelta(lt,'D'))
        obs = observations.sel(step=pd.Timedelta(lt,'D'))
        cm  = clim_mean.sel(step=pd.Timedelta(lt,'D'))

        subplots = fg(mod,dim)
        latex.set_style(style='white')
        fig,axes = plt.subplots(subplots[0],subplots[1],\
            figsize=latex.set_size(width='thesis',subplots=subplots,fraction=0.95),\
            subplot_kw=dict(projection=ccrs.PlateCarree()))

        axes = axes.flatten()

        x_group = list(mod.groupby(dim))
        y_group = list(obs.groupby(dim))
        cm_group = list(cm.groupby(dim))

        for n,(xlabel,xdata) in enumerate(x_group):

            ylabel,ydata   = y_group[n]
            cmlabel,cmdata = cm_group[n]

            # this approach is not bulletproof
            # check manually that right groups go together
            print('\tgraphics.skill_map: matched groups ',xlabel,ylabel)

            xdata  = xdata.unstack().sortby(['time'])
            ydata  = ydata.unstack().sortby(['time'])
            cmdata = cmdata.unstack().sortby(['time'])

            xdata,ydata,cmdata = xr.align(xdata,ydata,cmdata)

            try:
                ens_center_estimate = xdata.mean('member',skipna=True)
            except ValueError:
                ens_center_estimate = xdata

            try:
                ens_center_estimate_cm = cmdata.mean('member',skipna=True)
            except ValueError:
                ens_center_estimate_cm = cmdata

            score_mean   = xs.mae(ens_center_estimate,ydata,dim=[])
            score_clim   = xs.mae(ens_center_estimate_cm,ydata,dim=[])

            skill_score = ss(score_mean,score_clim,time='before',est='median')

            x = skill_score.lon.values
            y = skill_score.lat.values
            c = skill_score.values

            ax = axes[n]

            ax.coastlines(resolution='10m', color='black',\
                                    linewidth=0.2)

            ax.set_extent((0,25,55,75),crs=ccrs.PlateCarree())

            # cmap   = latex.cm_rgc(c='pink').reversed()
            cmap   = mpl.colors.ListedColormap(
                ['white','white','white','red','lightblue','royalblue','blue']
                                    )
            levels = [-1,-0.5,-0.25,-0.05,0.05,0.25,0.5,1]
            norm   = BoundaryNorm(levels,cmap.N)

            cs = ax.scatter(x, y, marker='o', c=c, s=1.2, cmap=cmap, norm=norm,
                alpha=0.8, transform=ccrs.PlateCarree())

            ax.set_title(month(xlabel))

        fig.colorbar(cs,ax=axes,boundaries=levels)
        fig.suptitle(title+' MAE at lead time: '+str(lt))
        save_fig(fig,fname)

def cluster_map(
              observations,
              model,
              clim_mean,
              clim_std,
              c_lim=(0.5,1),
              dim='validation_time.month',
              title='',
              lead_time=[9,16,23,30,37],
              filename=''
             ):

    for lt in lead_time:

        fname = 'mae_cluster_map/'+filename+'_'+dim.split('.')[1]+'_'+str(lt)

        mod = model.sel(step=pd.Timedelta(lt,'D'))
        obs = observations.sel(step=pd.Timedelta(lt,'D'))
        cm  = clim_mean.sel(step=pd.Timedelta(lt,'D'))
        cs  = clim_std.sel(step=pd.Timedelta(lt,'D'))

        subplots = fg(mod,dim)
        latex.set_style(style='white')
        fig,axes = plt.subplots(subplots[0],subplots[1],\
            figsize=latex.set_size(width='thesis',subplots=subplots,fraction=0.95),\
            subplot_kw=dict(projection=ccrs.PlateCarree()))

        axes = axes.flatten()

        x_group = list(mod.groupby(dim))
        y_group = list(obs.groupby(dim))
        cm_group = list(cm.groupby(dim))
        cs_group = list(cs.groupby(dim))

        for n,(xlabel,xdata) in enumerate(x_group):

            ylabel,ydata   = y_group[n]
            cmlabel,cmdata = cm_group[n]
            cslabel,csdata = cs_group[n]

            # this approach is not bulletproof
            # check manually that right groups go together
            print('\tgraphics.mae.cluster_map: matched groups ',xlabel,ylabel)

            xdata  = xdata.unstack().sortby(['time'])
            ydata  = ydata.unstack().sortby(['time'])
            cmdata = cmdata.unstack().sortby(['time'])
            csdata = csdata.unstack().sortby(['time'])

            xdata,ydata,cmdata,csdata = xr.align(xdata,ydata,cmdata,csdata)

            score_mean   = xs.mae(xdata.mean('member',skipna=True),ydata,dim=[])
            score_clim   = xs.mae(cmdata,ydata,dim=[])

            x,y,c = [],[],[]
            for loc in score_mean.location.values:

                target = score_mean.sel(location=loc)

                x.append(target.lon.values)
                y.append(target.lat.values)

                sm = lc.cluster(score_mean,loc,c_lim[0],c_lim[1])\
                                .mean('location',skipna=True)
                sc = lc.cluster(score_clim,loc,c_lim[0],c_lim[1])\
                                .mean('location',skipna=True)

                skill_score = ss(sm,sc,time='before',est='median')

                c.append(skill_score.values)

            ax = axes[n]

            ax.coastlines(resolution='10m', color='black',\
                                    linewidth=0.2)

            ax.set_extent((0,25,55,75),crs=ccrs.PlateCarree())

            # cmap   = latex.cm_rgc(c='pink').reversed()
            cmap   = mpl.colors.ListedColormap(
                ['white','white','white','red','lightblue','royalblue','blue']
                                    )
            levels = [-1,-0.5,-0.25,-0.05,0.05,0.25,0.5,1]
            norm   = BoundaryNorm(levels,cmap.N)

            cs = ax.scatter(x, y, marker='o', c=c, s=1.2, cmap=cmap, norm=norm,
                alpha=0.8, transform=ccrs.PlateCarree())

            ax.set_title(month(xlabel))

        fig.colorbar(cs,ax=axes,boundaries=levels)
        fig.suptitle(title+' MAE at lead time: '+str(lt))
        save_fig(fig,fname)
