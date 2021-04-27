from .handle_domain import update_storage_dict

domainID    = 'NVK'
model       = 'era'
ftype       = 'reanalysis'
var_abb     = 'sst'
t_res       = 'monthly'
i_time      = (1993,1) # start time, fmt: (year,month)
e_time      = (2020,9) # end time, fmt: (year,month)
dtype       = 'float64'
store_path  = ''

update_storage_dict(
                    domainID,
                    ftype,
                    model,
                    var_abb,
                    t_res,
                    i_time,
                    e_time,
                    dtype,
                    store_path
                    )
