
from S2S.data_handler import ERA5, ECMWF_S2SF, ECMWF_S2SH
from S2S.process import Hindcast

bounds = (5,10,60,65)
variables    = ['sst']

t_start  = (2020,4,10)
t_end    = (2020,4,20)

for var in variables:
    for high_res in [False]:
        for loader,label in zip([ECMWF_S2SH],['H']):

            print('\n')
            print('Variable:',var)
            print('High res:',high_res)
            print('Type:',label)

            try:

                Hindcast(
                        var,
                        t_start,
                        t_end,
                        bounds,
                        high_res=high_res,
                        process=False,
                        download=False,
                        split_work=False
                    )

                print('Success')

            except (KeyError,ValueError) as e:
                print(repr(e))
