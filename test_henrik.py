import pandas as pd

from S2S.data_handler import ERA5, BarentsWatch
from S2S.process import Hindcast, Observations, Grid2Point

domainID = 'norwegian_coast'
var      = 'sst'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)


high_res = False
steps    = pd.to_timedelta([9,16,23,30,37,44],'D')

hindcast = Hindcast(
                        var,
                        t_start,
                        t_end,
                        domainID,
                        high_res=high_res,
                        steps=steps,
                        process=False,
                        download=True
                    )

# observations must be weekly mean values with a time dimension
observations = BarentsWatch().load(
                                    [
                                        'Hisdalen',
                                        'Stokkvika'
                                    ]
                                )[var]

observations = Observations(
                            name='BarentsWatch',
                            observations=observations,
                            forecast=hindcast,
                            process=False
                            )

point_hindcast = Grid2Point(observations,hindcast).correlation(
                                                        step_dependent=True
                                                            )


# import S2S.graphics.graphics as gr
#
# gr.timeseries()

# # import scripts.Henrik.create_domain_file
# # import S2S.organize_IMR
#
# # from S2S.organize_barentzwatch import organize_files
# # organize_files()
#
# import scripts.Henrik.process_point_location
#
# # import scripts.Henrik.vpf_sb
# # import scripts.Henrik.score_pf
#
# # import scripts.Henrik.process_hindcast_wind
#
# # import scripts.Henrik.vpf_rawEC
# # import scripts.Henrik.vpf_ppEC_test
# # import scripts.Henrik.vpf_combo_ERA
#
# # import scripts.Henrik.aggregate_skill
