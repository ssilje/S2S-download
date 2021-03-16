import gridpp
import numpy as np

def make_grid_object(grb_data):
    latlats, lonlons = np.meshgrid(
        grb_data.latitude.data, grb_data.longitude.data
    )
    grid_object = gridpp.Grid(latlats, lonlons)

    return grid_object