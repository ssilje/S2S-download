import unittest

from S2S.date_helpers import get_forcast_date_cycle

class TestHelpers(unittest.TestCase):

    def test_datecycle(self):
        cycle = get_forcast_date_cycle(
            start_year=2019,
            start_month=7,
            start_day=1,
            num_periods=2
        )

        self.assertEqual(cycle[0].strftime('%Y-%m-%d'), '2019-07-01')
        self.assertEqual(cycle[1].strftime('%Y-%m-%d'), '2019-07-08')