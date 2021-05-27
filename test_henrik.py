# import scripts.Henrik.create_domain_file
# import scripts.Henrik.download_imr

# import scripts.Henrik.organize_barentzwatch as ob
# ob.organize_files()

# import scripts.Henrik.get_eide

import scripts.Henrik.verify_point_forecast_clean
import scripts.Henrik.map
# import scripts.Henrik.verify_point_forecast_onbw

# import scripts.Henrik.test_domain
# import scripts.Henrik.geographic_overview

# import scripts.Henrik.figure_receipts
# import scripts.Henrik.s2s2local

# import xarray as xr
# path = '/nird/projects/NS9853K/DATA//SFE/hindcast/ECMWF/sst/sst_CY46R1_2020-05-18_hindcast.grb'
# data = xr.open_dataset(local_path+filename[0], engine='cfgrib')
# print(data)

# from scripts.Henrik.data_handler import ERA5,ECMWF_S2SH
#
# domainID = 'NVK'
#
# t_start  = (2020,1,23)
# t_end    = (2021,1,14)
#
# clim_t_start  = (2000,1,1)
# clim_t_end    = (2004,1,4)
#
#
# for key in ERA5().load('sst',clim_t_start,clim_t_end,domainID).dims:
#     print(key)
#
# print(ECMWF_S2SH().load('sst',t_start,t_end,domainID))
