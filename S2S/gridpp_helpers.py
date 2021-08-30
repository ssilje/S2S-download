import gridpp
import numpy as np

def make_grid_from_grb(grb_data):
    latlats, lonlons = np.meshgrid(
        grb_data.latitude.data, grb_data.longitude.data
    )
    grid_object = gridpp.Grid(latlats, lonlons)

    return grid_object


def make_points_from_grb(grb_data, valid_points):

    latlats, lonlons = np.meshgrid(
        grb_data.latitude.data, grb_data.longitude.data
    )

    grid_object = gridpp.Points(latlats[valid_points], lonlons[valid_points])

    return grid_object