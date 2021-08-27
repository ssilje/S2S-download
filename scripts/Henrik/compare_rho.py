import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from S2S.local_configuration import config
from S2S.location_cluster import loc_from_name as lfn
from S2S.graphics import latex, graphics

path = config['VALID_DB']
l_pars = {'Hisdalen':'Hisdalen','Langøy':'Langøy S'}

for location in ['Hisdalen','Langøy']:

    data = xr.open_dataset(path+location+'_rho.nc').to_array().rename('rho')\
                                                    .squeeze().drop('variable')

    locs = data.location.values
    loc  = str(lfn(l_pars[location]))

    data = data.sel(location=loc)

    latex.set_style(style='white')
    fig,axes = plt.subplots(3,4,
                        figsize=latex.set_size(width='thesis',subplots=(1,1)))
    axes = axes.flatten()

    for month in np.arange(1,13,1):

        dm = data.where(data.time.dt.month==month)
        ax = axes[month-1]
        ax.plot(dm.step.dt.days.values,dm.mean('time').values)
        ax.set_title(graphics.month(month))
        ax.set_ylim(-1,1)

    fig.suptitle(location)
    graphics.save_fig(fig,location+'_rho')
