
from S2S.data_handler import ERA5, ECMWF_S2SF, ECMWF_S2SH
from S2S.process import Hindcast

bounds = (5,10,60,65)
variables    = ['u10']

t_start  = (2020,4,10)
t_end    = (2020,5,20)

for var in variables:
    for high_res in [False]:
        for loader,label in zip([ECMWF_S2SH],['H']):

            data = Hindcast(
                    var,
                    t_start,
                    t_end,
                    bounds,
                    high_res=high_res,
                    process=True,
                    download=False,
                    split_work=False,
                    cross_val=False
                )

            print(data.mean)
            print(data.std)
