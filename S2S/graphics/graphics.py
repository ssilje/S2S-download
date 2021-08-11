import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
import pandas as pd
import xskillscore as xs
import json
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm

from S2S.local_configuration import config
import S2S.xarray_helpers as xh

from . import latex

def quick_map(da,point=None):

    try:
        da = da.isel(member=0)
    except ValueError:
        pass
    try:
        da = da.isel(time=0)
    except ValueError:
        pass
    try:
        da = da.isel(step=0)
    except ValueError:
        pass

    p = da.transpose('lat','lon').plot(
            subplot_kws=dict(projection=ccrs.PlateCarree(),
            facecolor="white"),
            transform=ccrs.PlateCarree(),
        )
    p.axes.coastlines()

    if point:
        lat = point[1]
        lon = point[0]

        p.axes.scatter(lon, lat, marker='o',c='k',s=40,
            alpha=1, transform=ccrs.PlateCarree())
    plt.show()
    plt.close()

def loc_from_name(loc):
    """
    Returns name of Barentswatch location from loc number
    """
    try:

        with open(config['SITES'], 'r') as file:
            data = json.load(file)
            for line in data:
                if line['name']==loc:
                    return line["localityNo"]

    except ValueError:
        return loc

def name_from_loc(loc):
    """
    Returns name of Barentswatch location from loc number
    """
    try:
        _ = int(loc)

        with open(config['SITES'], 'r') as file:
            data = json.load(file)
            for line in data:
                if line["localityNo"]==int(loc):
                    return line['name']

    except ValueError:
        return loc


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
    if len(filename.split('/'))==2:
        path = filename.split('/')[0]
        make_dir(config['SAVEFIG']+'/'+path)
    else:
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
    try:
        ii = int(ii)
        ii = ['JAN','FEB','MAR',\
            'APR','MAI','JUN',\
            'JUL','AUG','SEP',\
            'OCT','NOV','DES'][ii-1]
    except (ValueError,IndexError) as error:
        pass
    return ii


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

    # if location dimension does not exist, assign it
    if np.isin(np.array(dax_in.dims),'location').sum() == 0:
        dax_in = dax_in.expand_dims('location')
        day_in = day_in.expand_dims('location')

    if dax_in.location.values.ndim==0:
        locations = [dax_in.location]

    else:
        locations = dax_in.location

    for loc in locations:
        for lt in dax_in.step:

            dax = dax_in.sel(location=loc,step=lt)
            day = day_in.sel(location=loc,step=lt)

            dax,day  = xr.broadcast(dax,day)

            if dim=='validation_time.month':
                dax = xh.assign_validation_time(dax)
                day = xh.assign_validation_time(day)

            fname    = 'qqplot_'+filename+'_'+name_from_loc(str(dax.location.values))+str(lt.dt.days.values)
            suptitle = title + ' ' + \
                            name_from_loc(str(dax.location.values)) +' lead time: '+str(lt.dt.days.values)#+\
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

            limits = [-5,5]

            x_group = list(dax.groupby(dim))
            y_group = list(day.groupby(dim))

            for n,(xlabel,xdata) in enumerate(x_group):

                ylabel,ydata = y_group[n]
                # this approach is not bulletproof
                # check manually that right groups go together
                print('\tgraphics.qq_plot: matched groups ',xlabel,ylabel)

                x,y = xr.align(xdata,ydata)
                x,y = x.values.flatten(),y.values.flatten()
                all = np.array([x,y]).flatten()

                ax = axes[n]

                ax.set_title(month(xlabel))

                ax.set_xlabel(x_axlabel)
                ax.set_ylabel(y_axlabel)

                ax.set_aspect('equal')
                ax.set_xlim(limits)
                ax.set_ylim(limits)
                ax.plot([limits[0], limits[1]], [limits[0], limits[1]],
                                                'k',alpha=0.7,linewidth=0.6)

                ax.plot(x,y,'o',alpha=0.4,ms=1)

            fig.suptitle(suptitle)
            save_fig(fig,fname)

def skill_plot(in_mod,in_clim,dim='validation_time.month',
                                filename='',
                                title='',
                                ylab=''):

    for loc in in_mod.location:

        mod  = in_mod.sel(location=loc)
        clim = in_clim.sel(location=loc)

        if dim.split('.')[0]=='validation_time':
            mod = xh.assign_validation_time(mod)
            clim = xh.assign_validation_time(clim)

        fname    = 'skill/'+filename+'_'+\
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
            print('\tgraphics.skill_plot: matched groups ',xlabel,ylabel)

            xdata = xdata.unstack().sortby(['time','step'])
            ydata = ydata.unstack().sortby(['time','step'])

            xdata,ydata = xr.align(xdata,ydata)

            ss = 1-(xdata/ydata).mean('time',skipna=True)

            zx = xdata.mean('time',skipna=True)
            zy = ydata.mean('time',skipna=True)

            x = np.array([td.days for td in ss.step.to_pandas()])
            z = ss.values

            ax = axes[n]

            ax.set_title(month(xlabel))
            ax.set_xlabel('lead time [D]')
            ax.set_ylabel(ylab)

            ax.plot(x,z,'o-',alpha=0.4,ms=1,label='ss')
            ax.plot(x,zx,'o-',alpha=0.4,ms=1,label='fc')
            ax.plot(x,zy,'o-',alpha=0.4,ms=1,label='clim')
            ax.legend()
            # ax.set_ylim((-1,1))
            ax.plot(
                    [0, 1],[0.5, 0.5],'k',
                    transform=ax.transAxes,alpha=0.7,
                    linewidth=0.6
                    )

        fig.suptitle(suptitle)
        save_fig(fig,fname)

def residual_plot(in_mod,in_clim,dim='validation_time.month',
                                filename='',
                                title='',
                                ylab=''):

    for loc in in_mod.location:

        mod  = in_mod.sel(location=loc)
        clim = in_clim.sel(location=loc)

        if dim.split('.')[0]=='validation_time':
            mod = xh.assign_validation_time(mod)
            clim = xh.assign_validation_time(clim)

        fname    = 'residuals/'+filename+'_'+\
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
            print('\tgraphics.skill_plot: matched groups ',xlabel,ylabel)

            xdata = xdata.unstack().sortby(['time','step'])
            ydata = ydata.unstack().sortby(['time','step'])

            xdata,ydata = xr.align(xdata,ydata)

            xx = xdata
            yy = ydata

            x = np.array([td.days for td in xx.step.to_pandas()])
            xx = xx.values.transpose()
            yy = yy.values.transpose()

            ax = axes[n]

            ax.set_title(month(xlabel))
            ax.set_xlabel('lead time [D]')
            ax.set_ylabel(ylab)

            ax.plot(x,z,'o',color='k',alpha=0.4,ms=1)
            # ax.plot(x,zx,'o-',alpha=0.4,ms=1,label='fc')
            # ax.plot(x,zy,'o-',alpha=0.4,ms=1,label='clim')
            # ax.legend()
            # ax.set_ylim((-1,1))
            ax.plot(
                    x,np.full_like(x,0),'k',
                    alpha=0.7,
                    linewidth=0.6
                    )

        fig.suptitle(suptitle)
        save_fig(fig,fname)

def timeseries(
            observations,
            cast,
            lead_time=[9,16,23],
            olab='obs',
            clabs=['forecast'],
            filename='',
            title=''
            ):
    """
    Plots timeseries of observations and forecast/hindcast

    args:
        cast list of xarray.DataArray
    """
    observations = xh.assign_validation_time(observations)

    if observations.location.values.ndim==0:
        locations = [observations.location]

    else:
        locations = observations.location

    for loc in locations:

        latex.set_style(style='white')
        fig,axes = plt.subplots(len(lead_time),1,
                        figsize=latex.set_size(width='thesis',
                            subplots=(1,1))
                        )
        clabs = clabs*len(cast)

        fname = 'timeseries/timeseries_'+\
                    filename+'_'+str(name_from_loc(loc.values))

        for n,lt in enumerate(lead_time):

            fname += '_'+str(lt)

            if hasattr(axes,'__iter__'):
                ax = axes[n]
            else:
                ax = axes


            ax.set_xlabel('Forecasted time')
            ax.set_ylabel('')

            try:
                o   = observations\
                        .sel(step=pd.Timedelta(lt,'D'),location=loc)\
                            .sortby('time')
            except ValueError:
                o   = observations\
                        .sel(step=pd.Timedelta(lt,'D'))\
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

                c = xh.assign_validation_time(c)

                try:
                    c = c\
                            .sel(step=pd.Timedelta(lt,'D'),location=loc)\
                                .sortby('time')
                except ValueError:
                    c = c\
                            .sel(step=pd.Timedelta(lt,'D'))\
                                .sortby('time')

                c = c.sel(time=slice(o.time.min(),o.time.max()))

                try:
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

                except ValueError:

                    ax.plot(
                            c.validation_time.squeeze(),
                            c.squeeze(),
                            alpha=0.7,
                            label=clabs[cn],
                            linewidth=0.5,
                            zorder=30
                        )

                    ax.fill_between(
                            [],
                            [],
                            []
                        )

                    score = xs.mae(o,c,dim=[])\
                                .mean('time',skipna=True)

                subtitle += ' '+clabs[cn]+' '+str(round(float(score.values),2))

            ax.set_title('Lead time '+str(lt)+' days '+subtitle)
            ax.legend()

        fig.suptitle(str(name_from_loc(loc.values))+' '+title)
        save_fig(fig,fname)

def point_map(da,c='r',poi=None,bw_special=False,boxes=None,bw_poi=None):

    name = poi

    latex.set_style(style='white')
    fig,axes = plt.subplots(1,1,\
        figsize=latex.set_size(width=345,subplots=(1,1),fraction=0.95),\
        subplot_kw=dict(projection=ccrs.PlateCarree()))

    ax = axes

    ax.coastlines(resolution='10m', color='grey',\
                            linewidth=0.2)

    ax.set_extent((0,25,55,75),crs=ccrs.PlateCarree())


    for loc in da.location:

        z = da.sel(location=loc)
        x = z.lon.values
        y = z.lat.values

        ax.plot(x, y, marker='o', color='orange', markersize=0.5,
            alpha=0.9, transform=ccrs.PlateCarree())
    if poi:
        poi = loc_from_name(poi)
        z = da.sel(location=str(poi))
        x = z.lon.values
        y = z.lat.values

        ax.plot(x, y, marker='o', color=c, markersize=1,
            alpha=1, transform=ccrs.PlateCarree())

        ax.text(x+1., y, name,
            verticalalignment='center', horizontalalignment='left',
            transform=ccrs.PlateCarree(),
            bbox=dict(facecolor='red', alpha=0.5, boxstyle='round'))

    if bw_special:

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.3, linestyle='--')

        for coords,label in boxes:
            box(ax,coords,proj=ccrs.PlateCarree(),label=label)

        if bw_poi is not None:
            for loc in bw_poi.location:

                z = da.sel(location=loc)
                x = z.lon.values
                y = z.lat.values

                ax.plot(x, y, marker='x', color='blue', markersize=0.3,
                    alpha=0.9, transform=ccrs.PlateCarree())

        ax.set_title('Fish farm locations')

        save_fig(fig,'fish_farm_loctions')

def box(ax,coords,proj,color='k',lwd=0.5,label=''):

    if len(coords)==4:
        x1 = coords[0]
        x2 = coords[1]
        y1 = coords[2]
        y2 = coords[3]

        ax.plot([x1,x1],[y1,y2],color=color,transform=proj,linewidth=lwd)
        ax.plot([x2,x2],[y1,y2],color=color,transform=proj,linewidth=lwd)
        ax.plot([x1,x2],[y1,y1],color=color,transform=proj,linewidth=lwd)
        ax.plot([x1,x2],[y2,y2],color=color,transform=proj,linewidth=lwd)

        ax.text(x1+0.1, y2, label,
            verticalalignment='top', horizontalalignment='left',
            transform=proj
            # ,
            # bbox=dict(facecolor='red', alpha=1, boxstyle='round')
            )
    if len(coords)==6:
        x1 = coords[0]
        x2 = coords[1]
        x3 = coords[2]
        x4 = coords[3]
        y1 = coords[4]
        y2 = coords[5]

        ax.plot([x1,x3],[y1,y2],color=color,transform=proj,linewidth=lwd)
        ax.plot([x2,x4],[y1,y2],color=color,transform=proj,linewidth=lwd)
        ax.plot([x1,x2],[y1,y1],color=color,transform=proj,linewidth=lwd)
        ax.plot([x3,x4],[y2,y2],color=color,transform=proj,linewidth=lwd)

        ax.text(x3+0.1, y2, label,
            verticalalignment='top', horizontalalignment='left',
            transform=proj
            # ,
            # bbox=dict(facecolor='red', alpha=1, boxstyle='round')
            )

    if len(coords)==8:
        x1 = coords[0]
        x2 = coords[1]
        x3 = coords[2]
        x4 = coords[3]
        y1 = coords[4]
        y2 = coords[5]
        y3 = coords[6]
        y4 = coords[7]

        ax.plot([x1,x3],[y1,y2],color=color,transform=proj,linewidth=lwd)
        ax.plot([x2,x4],[y3,y4],color=color,transform=proj,linewidth=lwd)
        ax.plot([x1,x2],[y1,y3],color=color,transform=proj,linewidth=lwd)
        ax.plot([x3,x4],[y2,y4],color=color,transform=proj,linewidth=lwd)

        ax.text(x3+0.1, y2, label,
            verticalalignment='top', horizontalalignment='left',
            transform=proj
            # ,
            # bbox=dict(facecolor='red', alpha=1, boxstyle='round')
            )

def skill_map(
              model,
              climate,
              dim='validation_time.month',
              title='SS',
              lead_time=[9,16,23],
              filename=''
             ):

    for lt in lead_time:

        fname = filename+'_'+str(lt)
        mod  = model.sel(step=pd.Timedelta(lt,'D'))
        clim = climate.sel(step=pd.Timedelta(lt,'D'))

        subplots = fg(mod,dim)
        latex.set_style(style='white')
        fig,axes = plt.subplots(subplots[0],subplots[1],\
            figsize=latex.set_size(width='thesis',subplots=subplots,fraction=0.95),\
            subplot_kw=dict(projection=ccrs.PlateCarree()))

        axes = axes.flatten()

        x_group = list(mod.groupby(dim))
        y_group = list(clim.groupby(dim))

        for n,(xlabel,xdata) in enumerate(x_group):

            ylabel,ydata = y_group[n]

            # this approach is not bulletproof
            # check manually that right groups go together
            print('\tgraphics.skill_map: matched groups ',xlabel,ylabel)

            xdata = xdata.unstack().sortby(['time'])
            ydata = ydata.unstack().sortby(['time'])

            xdata,ydata = xr.align(xdata,ydata)

            ss = 1-xdata.mean('time',skipna=True)/ydata.mean('time',skipna=True)

            x = ss.lon.values
            y = ss.lat.values
            c = ss.values

            ax = axes[n]

            ax.coastlines(resolution='10m', color='black',\
                                    linewidth=0.2)

            ax.set_extent((0,25,55,75),crs=ccrs.PlateCarree())

            cmap   = latex.cm_rgc(c='yellow').reversed()
            levels = np.arange(-1,1.2,0.2)
            norm   = BoundaryNorm(levels,cmap.N)

            cs = ax.scatter(x, y, marker='o', c=c, s=1, cmap=cmap, norm=norm,
                alpha=0.7, transform=ccrs.PlateCarree())

            ax.set_title(month(xlabel))

        fig.colorbar(cs,ax=axes)
        fig.suptitle(title+' for lead time: '+str(lt))
        save_fig(fig,fname)
