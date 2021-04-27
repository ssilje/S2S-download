import numpy as np

def arg_min_pos_val(array):
    """
    args:
        array: numpy array

    returns:
        index: int
    """
    return np.where(array >= 0., array, np.inf).argmin()

def arg_max_neg_val(array):
    """
    args:
        array: numpy array

    returns:
        index: int
    """
    return np.where(array <= 0., array, -np.inf).argmax()

def cap_grid(lon,lat,var,bounds):
    # TODO: make function able to interfere with boundaries of lon/lat
    """
    Caps grid at coords specified in bounds

    args:
        lon:    array of float (1-dim)
        lat:    array of float (1-dim)
        var:    array of float (3-dim)
        bounds: tuple of float (4-inst)

    returns:
        lon, lat, var: numpy.array (1,1,3-dim)

    - Format bounds: (min lon, max lon, min lat, max lat)

    - Includes the nearest grid point on the boundaries or outside bounds

    - Works with latitude in both increasing and decreasing order.
      Returns latitud in increasing order either way.
    """
    min_lon, max_lon, min_lat, max_lat = bounds

    lon = np.array(lon)
    lat = np.array(lat)
    var = np.array(var)

    min_lon_idx = arg_max_neg_val(lon-min_lon)
    max_lon_idx = arg_min_pos_val(lon-max_lon)

    min_lat_idx = arg_max_neg_val(lat-min_lat)
    max_lat_idx = arg_min_pos_val(lat-max_lat)

    # print(
    #     'Grid capped at',
    #     lon[min_lon_idx],
    #     lon[max_lon_idx],
    #     lat[min_lat_idx],
    #     lat[max_lat_idx],
    #     )

    # latitude organized in decreasin order, flip befor return
    if lat[0]>lat[-1]:

        return  lon[min_lon_idx:max_lon_idx+1],\
                np.flip(lat[max_lat_idx:min_lat_idx+1]),\
                np.flip(var[...,
                        max_lat_idx:min_lat_idx+1,
                        min_lon_idx:max_lon_idx+1],
                    axis=-2)

    # latitude organzed in increasing order, return current
    else:

        return  lon[min_lon_idx:max_lon_idx+1],\
                lat[min_lat_idx:max_lat_idx+1],\
                var[...,min_lat_idx:max_lat_idx+1,min_lon_idx:max_lon_idx+1]
