"""Paths and methods for boundary files (uas, csa, counties, blocks, etc)"""
from utilities import measurements

csa_path = "data/cb_2017_us_csa_500k/cb_2017_us_csa_500k.shp"
cbsa_path =  "data/cb_2017_us_cbsa_500k/cb_2017_us_cbsa_500k.shp"
ua_path = "data/cb_2017_us_ua10_500k/cb_2017_us_ua10_500k.shp"

zip_boundaries_path = "data/cb_2017_us_zcta510_500k/cb_2017_us_zcta510_500k.shp"
county_boundaries_path = "data/cb_2017_us_county_500k/cb_2017_us_county_500k.shp"
texas_blocks_path = "data/ACS_2016_5YR_BG_48_TEXAS.gdb"
texas_places_path = measurements.texas_places_path