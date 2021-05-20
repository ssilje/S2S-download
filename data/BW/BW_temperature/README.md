# BarentzWatch export

Data delivered by BarentsWatch (https://www.barentswatch.no/fiskehelse/). Processed by Benjamin Narum.

Information from BarentsWatch APIs are made available by the Norwegian Licence for Open Government Data (NLOD) (https://data.norge.no/nlod/en) unless API documentation specify otherwise. See details on terms at https://www.barentswatch.no/en/about/open-data-via-barentswatch/.

## Description

BW_temperature_export.csv

- week
- year
- date = date of Wednesday representing a week
- locNo = Official site ID of locality
- temperature = reported temperature in Celsius
- active_status = whether the locality was active in a given week (temperature reports not are always given otherwise)

metadata_BW_sites.json = mapping between coordinates and locNo.

