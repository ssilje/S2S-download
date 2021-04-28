# environment dependencies
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import cartopy.crs as ccrs

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
    data.load('ERA5',(2020,1,1),(2020,1,2),domainID+'_map')

    sst = data.ERA5.sst
    lon = data.ERA5.lon
    lat = data.ERA5.lat

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(color='black', linewidth=0.2)
    # ax.set_extent((map_bounds),crs=ccrs.PlateCarree())
    cs = ax.contourf(lon,lat,sst.mean(axis=0),projection=ccrs.PlateCarree())
    mark_box(ax,domain_bounds,ccrs.PlateCarree(),lwd = 1.)
    plt.colorbar(cs)

    save_fig(domainID+'_map.png')

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
