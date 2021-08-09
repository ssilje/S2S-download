
from datetime import datetime

# dictionary of model versions with (start,end)
model_version_specs = dict(
    ECMWF = dict(
        CY43R1 = ('2016-11-22','2017-07-10'),
        CY43R3 = ('2017-07-11','2018-06-05'),
        CY45R1 = ('2018-06-06','2019-06-10'),
        CY46R1 = ('2019-06-11','2020-06-29'),
        CY47R1 = ('2020-06-30','2021-05-10'),
        CY47R2 = ('2021-05-11',datetime.strftime(datetime.today(),"%Y-%m-%d"))
    )
)


def which_mv_for_init(fc_init_date,model='ECMWF',fmt='%Y-%m-%d'):
    """
    return model version for a specified initialization date and model

    INPUT:
            fc_init_date:   date string YYYY-mm-dd
            model:          string for the modeling center (currently just
                            'ECMWF' is valid)
            fmt:            string specifying the date format,
                            default: '%Y-%m-%d'
    OUTPUT:
            model version as string
    """

    # convert date string to datetime object:
    fc_init_datetime = datetime.strptime(fc_init_date,fmt)

    # got through the model versions from the above dictionary:
    for MV,mv_dates in model_version_specs[model].items():
        # convert first and last dates to datetime:
        mv_first = datetime.strptime(mv_dates[0],fmt),
        mv_last  = datetime.strptime(mv_dates[-1],fmt)
        # check if the given date is within the current model version's
        # start and end dates:
        if  mv_first <= fc_init_datetime <= mv_last:
            valid_version = MV

    try:
        return valid_version
    except:
        print('No matching model version found...')
        return None
