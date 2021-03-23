import pandas as pd

def get_forcast_date_cycle(
        start_year,
        start_month,
        start_day,
        num_weeks,
):

    start_date = pd.to_datetime(f"{start_year:04}-{start_month:02}-{start_day:02}")

    if start_date.weekday() != 0:
        raise AttributeError("Date is not a Monday. Forecasts start Mondays.")

    dr_mondays = pd.date_range(
        start=start_date,
        periods=num_weeks,
        freq='W-MON'
    )
    return dr_mondays

    # dr_thursdays = pd.date_range(
    #     start=start_date,
    #     periods=num_weeks,
    #     freq='W-THU'
    # )

    # return dr_mondays.union(dr_thursdays)
