import numpy as np

# def distance(lat0,lat1,lon0,lon1):
#
#     R = 6378.137
#
#     dlat = lat1 * np.pi / 180 - lat0 * np.pi / 180
#     dlon = lon1 * np.pi / 180 - lon0 * np.pi / 180
#
#     a = np.sin(dlat/2) * np.sin(dlat/2) +\
#         np.cos(lat0 * np.pi / 180) * np.cos(lat1 * np.pi / 180) *\
#         np.sin(dlon/2) * np.sin(dlon/2)
#
#     c = 2 * np.atan2(np.sqrt(a), np.sqrt(1-a))

latitude     = 62.5
r            = 6378137
new_r        = r * np.sin( np.pi * latitude / 180 )

dist_per_rad = new_r

rad800 = 2 * np.pi * 800/new_r

xend   = 38.5
xstart = 0

yend   = 75
ystart = 50

inc = rad800/np.pi * 180
xsize = (xend-xstart)/inc
ysize = (yend-ystart)/inc

print('inc',inc)
print('xsize',xsize)
print('ysize',ysize)
