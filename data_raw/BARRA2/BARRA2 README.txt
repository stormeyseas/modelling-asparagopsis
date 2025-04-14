################################################################################

Bureau of Meteorology Atmospheric high-resolution Regional Reanalysis for 
Australia, BARRA Version 2

################################################################################

VERSION
v1.0: 20/12/23: Created Dr Chun-Hsu Su, chunhsu.su@bom.gov.au
v1.1: 27/08/24: Updated for AUST-04/BARRA-C2 and BARRA-R2 convective parameters
                data sets.

--------------------------------------------------------------------------------

INTRODUCTION

BARRA2 provides the Bureau's higher resolution regional atmospheric reanalysis 
over Australia and surrounding regions, spanning 1979-present day time period. 
When completed, it replaces the first version of BARRA (Su et al., 
doi: 10.5194/gmd-14-4357-2021; 10.5194/gmd-12-2049-2019).

It is produced using the Bureau's data assimilation system for numerical weather
prediction - 4D variational scheme, and ACCESS as a limited-area dynamical
coupled atmosphere-land model - Unified Model (UM) and JULES.

The data set includes sub-daily, daily and monthly data for temperature,
moisture, wind and flux variables at sub-surface, surface, and pressure levels,
and heights above surface. The vertical levels include many pressure levels and
several heights above surface.

Data Provider: Bureau of Meteorology

NCI Data Catalogue: https://dx.doi.org/10.25914/1x6g-2v48

NCI THREDDS Data Server: https://dx.doi.org/10.25914/1x6g-2v48

License: https://creativecommons.org/licenses/by/4.0/

Extended Documentation: https://opus.nci.org.au/x/DgDADw 

--------------------------------------------------------------------------------

CONDITIONS OF USE

The data collection is considered a research product containing direct modelling 
outputs from the ACCESS (Australian Community Climate and Earth system Simulator) 
that has not been fully evaluated. 

Users are advised to refer to the Extended Documentation for the known issues 
and undertake evaluation of the data prior to use in their applications. 

The Bureau of Meteorology seeks user feedback on the quality and usage of the 
data, to help identify areas for improvements. 

The Bureau can also advise appropriate use of the data. 

Please refer feedback and questions on data use to help@nci.org.au.

--------------------------------------------------------------------------------

DISCLAIMER

Refer to http://www.bom.gov.au/other/disclaimer.shtml

--------------------------------------------------------------------------------

CITATION

   Please cite both,
   Bureau of Meteorology. (2023). Bureau of Meteorology Atmospheric high-resolution Regional Reanalysis for Australia - Version 2 (BARRA2) (Version <INSERT>). NCI Australia. https://doi.org/10.25914/1X6G-2V48

   and

   (for BARRA-R2/RE2)
   Su, C.-H., Rennie, S., Dharssi, I., Torrance, J., Smith, A., Le, T., Steinle, P., Stassen, C., Warren, R. A., Wang, C., and Le Marshall, J. (2022), BARRA2: Development of the next-generation Australian regional atmospheric reanalysis, Bureau Research Report 067, accessed online http://www.bom.gov.au/research/publications/researchreports/BRR-067.pdf

   (for BARRA-C2)
   Su, C.-H., Rennie, S., Torrance, J., Howard, E., Stassen, C., Lipson, M., Warren, R., Pepler, A., Dharssi, I., Franklin, C. (2024), BARRA-C2: Development of the kilometre-scale downscaled atmospheric reanalysis over Australia, Bureau Research Report 097, accessed online http://www.bom.gov.au/research/publications/researchreports/BRR-097.pdf
   
--------------------------------------------------------------------------------

FILE ORGANISATION

   /g/data/ob53
   |-- <product>
     |-- <nature of data>
             |-- <activity_id>
                  |-- <domain_id>
                       |-- <RCM-institution_id>
                            |-- <driving_source_id>
                                 |-- <driving_experiment_id>
                                      |-- <driving_variant_label>
                                           |-- <source_id>
                                                |-- <version_realisation>
                                                     |-- <freq>
                                                          |-- <variable_id>
                                                              |-- <version>
                                                                   |-- <netcdf filename>

   where,
     <product> is BARRA2
     <nature of data> is output, referring to model output
     <activity_id> is reanalysis
     <domain_id> is spatial domain and grid resolution of the data, namely 
               AUS-11, AUS-22, AUST-04, AUST-11
     <RCM-institution> is BOM
     <driving_source_id> is ERA5 that drives BARRA2 at the lateral boundary
     <driving_experiment_id> is historical
     <driving_variant_label> labels the nature of ERA5 data used, either hres 
               or eda
     <source_id> is BARRA-R2, BARRA-RE2, or BARRA-C2, refer to 
               Extended Documentation
     <version_realisation> identifies the modelling version of BARRA2 (TBC on 
               identifying data version)
     <freq> is the time frequency of the data: 1hr (1-hourly), 3hr, 6hr, 
               day (daily), mon (monthly), fx (constant)
     <variable_id> is the variable name, mostly based on,
               https://docs.google.com/spreadsheets/d/1qUauozwXkq7r1g-L4ALMIkCNINIhhCPx/edit#gid=1672965248
     <version> denotes the date of data generation or date of data release
               or 'latest' points to the latest version.
     <netcdf filename> is
               <variable_id>_<domain_id>_<driving_source_id>_<driving_experiment_id>_<driving_variant_label>_<RCM-institution_id>_<source_id>_<version_realisation>_<freq>[_<StartTime>-<EndTime>].nc

  Example:
    /g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/mon/ua100m/v20231001/ua100m_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_mon_197901-197901.nc
    /g/data/ob53/BARRA2/output/reanalysis/AUS-22/BOM/ERA5/historical/eda/BARRA-RE2/v1/1hr/tas/v20231001/tas_AUS-22_ERA5_historical_eda_BOM_BARRA-RE2_v1_1hr_202201-202201.nc

--------------------------------------------------------------------------------

DATA FORMAT
netCDF-4 classic

--------------------------------------------------------------------------------

ACKNOWLEDGEMENTS

The production of BARRA2 was supported with funding from the Australian Climate 
Service. The availability of BARRA2 data at NCI has been made possible through 
the Australian Climate Service, and support of NCI.
