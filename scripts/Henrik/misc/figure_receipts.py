# environment dependencies
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import cartopy.crs as ccrs
import xarray as xr

# local dependencies
from S2S.local_configuration import config
from .create_domain_file import make_dir
from .handle_domain import get_bounds
from .data_handler import SST

def mark_box(ax,bounds,proj,col='gray',lwd = 0.75):
    """
    Mark subplot with box around bounds
    """
    lines = bounds
    ax.plot([lines[0],lines[0]],[lines[2],lines[3]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)
    ax.plot([lines[1],lines[1]],[lines[2],lines[3]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)
    ax.plot([lines[0],lines[1]],[lines[2],lines[2]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)
    ax.plot([lines[0],lines[1]],[lines[3],lines[3]],\
                    '-',color=col,transform=proj,linewidth=lwd,zorder=30)

def save_fig(filename):

    make_dir(config['SAVEFIG'])
    plt.savefig(config['SAVEFIG']+filename)
    plt.close()
    print('Figure stored at: '+config['SAVEFIG']+filename)

def geographic(domainID):

    domain_bounds = get_bounds(domainID)

    data = SST()
    data.load('ERA5',(2020,1,1),(2020,2,1),domainID+'_map')

    sst = np.array(data.ERA5.sst)
    lon = data.ERA5.lon
    lat = data.ERA5.lat

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(color='black', linewidth=0.2)
    # ax.set_extent((map_bounds),crs=ccrs.PlateCarree())
    cs = ax.contourf(lon,lat,sst.mean(axis=0),projection=ccrs.PlateCarree())
    mark_box(ax,domain_bounds,ccrs.PlateCarree(),lwd = 1.)
    plt.colorbar(cs)

    save_fig(domainID+'_map.png')

def reliability_plot(rel,sea,stp):
    # sns.set_theme(style="ticks")
    # fig, axes = plt.subplots(2,2)
    # for n,ax in enumerate(axes.flatten()):
    #
    #     for l in [10,20]:
    #
    #         data = rel[l][n].isel(threshold=0)
    #         ax.plot(data.forecast_probability,data.forecast_probability,'k')
    #         ax.plot(data.forecast_probability,data.interpolate_na(dim='forecast_probability').variable.squeeze(),label='pc: 0.75 lt:'+str(stp[l]))
    #
    #     # for l in [10,20]:
    #     #
    #     #     data = rel[l][n].isel(threshold=1)
    #     #     ax.plot(data.forecast_probability,data.forecast_probability,'k')
    #     #     ax.plot(data.forecast_probability,data.interpolate_na(dim='forecast_probability').variable.squeeze(),label='pc: 0.25 lt:'+str(stp[l]))
    #
    #     o = (data.variable * data.samples).sum(skipna=True)/data.samples.sum(skipna=True)
    #     o = np.float(o)
    #     ax.plot([0,1],[o,o],'--',color='black',label='no res')
    #
    #     fc_pb    = np.arange(0,1.1,0.001)
    #     no_skill_l = np.empty_like(fc_pb)
    #     no_skill_u = np.empty_like(fc_pb)
    #
    #     no_skill_u[(fc_pb+o)/2<o]  = (fc_pb+o)[(fc_pb+o)/2<o]/2
    #     no_skill_u[(fc_pb+o)/2>=o] = 1.
    #
    #     no_skill_l[(fc_pb+o)/2<o]  = 0
    #     no_skill_l[(fc_pb+o)/2>=o] = (fc_pb+o)[(fc_pb+o)/2>=o]/2
    #
    #     ax.fill_between(fc_pb, no_skill_l, no_skill_u,color='grey',alpha=0.2)
    #     # ax.plot(fc_pb,o+(fc_pb+0)/2,'-.',color='grey',label='no skill')
    #     ax.set_title(sea[n])
    #     ax.legend()
    #     ax.set_xlabel('P( y )')
    #     ax.set_ylabel('P( o | P(y) )')
    # plt.show()
    # plt.close()

    sns.set_theme(style="ticks")
    fig, axes = plt.subplots(2,2)
    for n,ax in enumerate(axes.flatten()):

        for l in [10,20]:

            data = rel[l][n].isel(threshold=0)
            x = np.transpose(np.array([data.forecast_probability.squeeze()]*2))
            y = np.array(data.samples.squeeze())
            y = np.transpose(np.array([np.zeros_like(y),y]))
            for a in range(len(x)):
                ax.plot(x[a]-(l/10-1.99),y[a],color=['blue','orange'][int(l/10)-1])
            ax.plot(x[a]-(l/10-1.99),y[a],label='pc: 0.75 lt:'+str(stp[l]),color=['blue','orange'][int(l/10-1)])

        ax.set_title(sea[n])
        ax.legend()
        ax.set_xlabel('P( y )')
    plt.show()
    plt.close()

def qq_plot(df,figure_name):

    sns.set_theme(style="ticks")
    sns.lmplot(x="observations", y="model", data=df,
                height=2, aspect=1.5, col='season',row='lead_time',
                scatter_kws={"s": 50, "alpha": 0.3})
    save_fig(figure_name+'_qq.png')

def line_plot(df,figure_name,var_name='values'):

    sns.set_theme(style="ticks")
    sns.relplot(
        data=df,
        x="lead_time", y=var_name, col="month",
        kind="line", palette="crest", linewidth=4, zorder=5,
        col_wrap=3, height=2, aspect=1.5)
    save_fig(figure_name+'_line.png')
