from S2S.data_handler import ERA5, ECMWF_S2SF, ECMWF_S2SH

bounds = (5,10,60,65)
variables    = ['sst','t2m','u10','v10']

t_start  = (2020,4,10)
t_end    = (2020,4,20)

for var in variables:
    for high_res in [True,False]:
        for loader,label in zip([ERA5,ECMWF_S2SF,ECMWF_S2SH],['ERA','H','F']):

            print('\n')
            print('Variable:',var)
            print('High res:',high_res)
            print('Type:',label)

            try:
                data = loader(high_res=high_res)\
                                .load(
                                        var,
                                        t_start,
                                        t_end,
                                        bounds
                                    )[var]

                print('Success')

            except (KeyError,ValueError) as e:
                print(repr(e))
