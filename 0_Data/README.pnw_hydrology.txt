----
<GCM>_<SCEN> and GMET_hist directories include:

- <GCM>_<SCEN>_mizuRoute_daily.nc: 
  desc:  
    mizuRoute daily mean flow at all the reaches in the river network (daily mean of 3hr raw output from mizuRoute)
  variables:
   local_runoff, (time, seg), m3/s, float, "local, lateral flow - remapped runoff times river model catchment area"
   streamflow,   (time, seg), m3/s, float, "total, routed flow using impulse response function reach routing method"
   reachID,      (seg),       -,    int64, "ID of river reach of MERIT-basin"

- <GCM>_<SCEN>_mizuRoute_daily_site.nc: 
  desc:  
    mizuRoute daily mean flow at 414 sites only. 
  variables:
   streamflow, (time, seg), m3/s, float,  "total, routed flow using impulse response function reach routing method"
   reachID,    (seg),       -,    int64,  "ID of river reach"
   site,       (seg),       -,    string, "flow site ID"

- <GCM>_<SCEN>_bmorph_site_univariate_daily.nc: 
  desc:  
    bmorph at-site bias corrected flow (independent bias correction, spatially consistent bias correction) and 
    raw simulated flow at 200 sites that were used for bmorph bias correction
  variables:
    flow_scbc_u, (time, site), m3/s", float,  "streamflow - spatially consistent bias correction"
   flow_ibc_u,   (time, site), m3/s,  float,  "streamflow - independent bias correction at flow site"
   flow_raw,     (time, site), m3/s,  float,  "streamflow not bias corrected"
   site,         (site),       -,     string, "flow site ID"

- <GCM>_<SCEN>_summa_daily.nc: 
  desc:
    summa daily outputs including swe, soil moisture, aquifer storage, canopy water, runoff, and ET. 
  variables:
    scalarTotalRunoff_mean, (time, hru), mm/day, float, "total runoff (mean)"
    scalarTotalET_mean,     (time, hru), mm/day, float, "total ET (mean)"
    scalarSWE,              (time, hru), mm,     float, "snow water equivalent (instant)"
    scalarCanopyWat,        (time, hru), mm,     float, "total water on the vegetation canopy (instant)"
    scalarTotalSoilWat,     (time, hru), mm,     float, "water in the soil (instant)"
    scalarAquiferStorage,   (time, hru), mm,     float, "storage of water in the aquifer (instant)"
    hru,                    (hru),       -,      int64, "ID defining the hydrologic response unit" 

- <GCM>_<SCEN>_summa_daily_basin_mean.nc: 
  desc:
    Daily mean summa variables averaged over upstream area for 414 sites.
  variables:
    scalarTotalRunoff_mean, (time, site), mm/day, float, "total runoff (mean)"
    scalarTotalET_mean,     (time, site), mm/day, float, "total ET (mean)"
    scalarSWE,              (time, site), mm,     float, "snow water equivalent (instant)"
    scalarTotalSoilWat,     (time, site), mm,     float, "water in the soil (instant)"
    scalarAquiferStorage,   (time, site), mm,     float, "storage of water in the aquifer (instant)"
    scalarCanopyWat,        (time, site), mm,     float, "total water on the vegetation canopy (instant)"
    site,                   (site),       -,      string, "naturalized site ID"

- <GCM>_<SCEN>_daily_t_p.nc: 
  desc:
    Daily Tmax, Tmin, and precipitation at each SUMMA HUC12.
  variables:
	  prec,  (time, hru), mm,       float, "daily precipitation"
	  t_min, (time, hru), C-degree, float, "daily minimum temperature"
	  t_max, (time, hru), C-degree, float, "daily maximum temperature"
	  hru,   (hru),       -,        int64, "hru ID"

- <GCM>_<SCEN>_daily_t_p_basin_mean.nc: 
  desc:
    Daily Tmax, Tmin, and precipitation averaged over upstream area for 414 sites. 
  variables:
	  prec,  (time, site), mm,       float,  "daily precipitation"
	  t_max, (time, site), C-degree, float,  "daily maximum temperature"
	  t_min, (time, site), C-degree, float,  "daily minimum temperature"
	  site,  (site),       -,        string, "naturalized site ID" 

- ancillary_data:
  geospatial_data: geopackages for MERIT-basin flowline/catchments, SUMMA HUC12 catchments and Naturalized flow sites. See ancillary_data/geospatial_data/NOTE for more details 
  PNW_unimpaired_flow_1951-2018.nc: naturalized flow data at 231 sites (some of sites do not have valid data, analysis used 214 sites. removed in PNW_unimpaired_flow_site.gpkg indicates which ones are used) 
