import numpy as np
import pandas as pd

def to_datetime(tup):
    """
    Returns pandas.DatetimeIndex from tuple

    args:
        tup: tuple of int, fmt: (year,month,date)

    returns:
        pandas.DatetimeIndex
    """
    return pd.to_datetime(f'{tup[0]}-{tup[1]}-{tup[2]}',format='%Y-%m-%d')

def forecast_cycle(t_start,t_end):
    """
    Returns pandas.DatetimeIndex with dates of a 7-day period and
    a 7 day period with a 4-day shift.

    args:
        t_start: tuple of int, fmt: (year,month,date)
        t_end:   tuple of int, fmt: (year,month,date)

    returns:
        pandas.DatetimeIndex
    """
    dates_fcycle_1 = pd.date_range(
                            start=f'{t_start[0]}-{t_start[1]}-{t_start[2]}',
                            end=f'{t_end[0]}-{t_end[1]}-{t_end[2]}',
                            freq='7D'
                            )

    dates_fcycle_2 = pd.date_range(
                            start=f'{t_start[0]}-{t_start[1]}-{t_start[2]+4}',
                            end=f'{t_end[0]}-{t_end[1]}-{t_end[2]+4}',
                            freq='7D'
                            )
    return pd.DatetimeIndex(\
                np.append(\
                    np.array(dates_fcycle_1),\
                    np.array(dates_fcycle_2)\
                    )\
                ).sort_values()


def weekly_forecast_cycle(t_start,t_end):
    """
    Returns pandas.DatetimeIndex with dates of a 7-day period.

    args:
        t_start: tuple of int, fmt: (year,month,date)
        t_end:   tuple of int, fmt: (year,month,date)

    returns:
        pandas.DatetimeIndex
    """
    return pd.date_range(
                    start=f'{t_start[0]}-{t_start[1]}-{t_start[2]}',
                    end=f'{t_end[0]}-{t_end[1]}-{t_end[2]}',
                    freq='7D'
                )


def days_from(t_start,t_end,match_fc=False):
    """
    Returns pandas.DatetimeIndex with dates of a 1-day period.

    args:
        t_start: tuple of int, fmt: (year,month,date)
        t_end:   tuple of int, fmt: (year,month,date)

    returns:
        pandas.DatetimeIndex
    """
    if not match_fc:
        return pd.date_range(
                    start=f'{t_start[0]}-{t_start[1]}-{t_start[2]}',
                    end=f'{t_end[0]}-{t_end[1]}-{t_end[2]}',
                    freq='D'
                    )
    else:
        start = to_datetime(t_start)
        end   = to_datetime(t_end) + pd.Timedelta(value=6,unit='W')
        return pd.date_range(start=start,end=end,freq='D')
