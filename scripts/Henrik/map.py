import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from scripts.Henrik.data_handler import BarentsWatch, ERA5, ECMWF_S2SH, Archive
import scripts.Henrik.xarray_helpers as xh
from S2S.local_configuration import config
import scripts.Henrik.models as models
import scripts.Henrik.organize_barentzwatch as ob
import scripts.Henrik.graphics as gr
import scripts.Henrik.latex as latex
from .handle_domain import get_bounds

import cartopy.crs as ccrs

domainID = 'full_grid'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2020,2,18)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

point_observations = BarentsWatch().load(['Hisdalen','Langøy S']).sortby('time')
observations = ERA5().load(var,clim_t_start,clim_t_end,domainID)[var]-272.15

domain_bounds = get_bounds(domainID)

sst = observations.mean('time').transpose('lat','lon')
lon = observations.lon
lat = observations.lat

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(color='black', linewidth=0.2)
cs = ax.contourf(lon,lat,sst,projection=ccrs.PlateCarree())
plt.colorbar(cs)
fig = plt.gca()
gr.save_fig(fig,'map')
