# S2S-download


------------------------------------------------------

### S2S: sub-seasonal to seasonal prediction project. 

Web: https://confluence.ecmwf.int/display/S2S

It is a WWRP/THORPEX-WCRP joint research project established to improve forecast skill and understanding on the sub-seasonal to seasonal time scale

The S2S data are downloaded from ECMWF's S2S database on the MARS archive: https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets

All the models are stored on the S2S-archive on a common 1.5/1.5 regular lat-lon grid (ocean parameters on 1x1 degree grid), and this is the resolution the data is downloaded as. 
Note that the different models have different original resolutions (and the model have different resolution for different lead-times).
See the S2S-documentation for more information: https://confluence.ecmwf.int/display/S2S/Models

The data is donwloaded globally as hindcast and forecast, with the associated ensemble members. 

------------------------------------------------------

The data structure: 
S2S/hindcast/_modeling-center_/_type-files_/_variable-name_/

where (add the information when new models, variables are included)
- _modeling-center_ : ECMWF
- _type-files_ : sfc, pl
- _variable-name_ : 
  - tp (accumulated): 228228
  - t2m: 167
  - sst: 34
  - mslp: 151
  
------------------------------------------------------
  
Information about the specific models that is downloaded: 

------------------------------------------------------
ECMWF: 
- Model version CY46R1, https://confluence.ecmwf.int/display/S2S/ECMWF+Model+Description+CY46R1)
- Data time of first forecast run:   11 June 2019
- First forecast downloaed: 2019-07-01 and one year (2020-06-29)
- Hindcast: 20 years
- Ensemble members: 
  - forecast 50 + 1
  - hindcast 10 + 1
