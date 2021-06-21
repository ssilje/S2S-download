import pandas as pd

from S2S.date_helpers import get_forcast_date_cycle

cycle = get_forcast_date_cycle(
    start_year=2019,
    start_month=7,
    start_day=1,
    num_periods=3
)