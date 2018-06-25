# Compiled by Asher Spector - these are not officially a part of the shapefiles

# Some notes on methodology, for consistency. I classify:
# (a) Areas which allow both multifamily and single family 
# dwelling as multifamily.
# (b) Areas which allow mixed-use development (i.e. commercial + residential)
# as some type of residential area in most cases. 
# (c) If multifamily housing is permitted but with some restrictions, 
# the area is classified as multifamily housing. 
# However, it is worth noting that that is not necessarily the case for the larger datasets which have their own methodologies
# (i.e. the Dallas + surrounding counties dataset). 

library(rgdal)
library(maptools)

# Inputs list -----------------------------------------------------------
# Name is the name of the municipality
# Path is the path to the shapefile
# sf is single family codes
# mf is multifamily codes
# or stands for other residential codes
# separator is the separator (str) by which to split codes into their component parts
# feature is the feature which actually lists the zoning codes per area
# -----------------------------------------------------------------------

# Basic projection inputs ----------------------------------------------- 

# For NAD83 - NO (HARN) - in ftUS
texas_north_central <- CRS('+init=EPSG:2276') #Dallas
texas_central <- CRS('+init=EPSG:2277') #Austin
texas_south_central <- CRS('+init=EPSG:2278') #Houston
lat_long <- CRS('+init=EPSG:4326') # standard lat and long


# AUSTIN -----------------------------------------------------------------

austin_inputs <- list(
  name = 'Austin',
  path = "austin_zoning/geo_export_8d172e2c-89a0-4c2f-a91d-337da740f107.shp",
  sf = 'SF', # The codes which correspond to single family
  mf = 'MF', # Codes which correspond to multifamily
  or = c('MH', 'RR', 'LA'), # Codes which correspond to other residential
  separator = '-',
  feature = 'zoning_zty',
  proj4string = lat_long
)

round_rock_inputs <- list(
  name = 'Round Rock', 
  path = "round_rock_zoning/ZONING_1.shp",
  feature = 'BASE_ZONIN',
  sf = c('SF', 'TH', 'SF1', 'SF2', 'SF3', 'SF4', 'SF5', 'SF6'),
  mf = c('MF', 'MU', 'TF'),
  or = c(),
  separator = '-',
  proj4string = texas_central
)

pflugerville_inputs <- list(
  name = 'Pflugerville',
  path = "pflugerville_zoning/Zoning_Districts.shp",
  feature = 'ZOINING_TY', # This is not a typo - it is spelled wrong in the data
  sf = c('A', 'SF'), 
  mf = c('2', 'CL3', 'CL4', 'CL5'),
  or= 'MH',
  separator = '-',
  proj4string = texas_central 
)

georgetown_inputs <- list(
  name = 'Georgetown',
  path = "georgetown_zoning/Zoning.shp",
  feature = 'ZONE', 
  sf = c('AG', 'RE', 'RL', 'RS', 'TH'),
  mf = c('TF', 'MF', 'MU', 'MU-DT', 'MUDT'),
  or = c('MH'),
  separator = '-',
  proj4string = texas_central
)

cedar_park_inputs <- list(
  name = 'Cedar Park',
  path = "cedar_park_zoning/Zoning__Zoning_Districts.shp",
  feature = 'ZoningType',
  separator = ' -',
  sf = c('RA', 'SR', 'SU', 'UR'),
  mf = c('MU', 'MF'),
  or = c(),
  proj4string = texas_central
)

hutto_inputs <- list(
  name = 'Hutto',
  path = "hutto_zoning/Zoning_Districts.shp",
  feature = 'ZONING',
  separator = 'no_separator',
  # Not totally sure about urban residential/residential classifications
  sf = c('Single Family', 'Residential'),
  mf = c('Two Family', 'Multi-Family', 'Urban Residential', 'Co-op District'),
  or = c(),
  proj4string = texas_central
)

# Need a data dictionary for this, otherwise it's useless
all_austin_inputs <- list(
  name <- 'Austin Area',
  path = 'central_land_use_2010/land_use_2010.shp',
  proj4string = texas_central
)

# DALLAS -----------------------------------------------------------------


# Downloaded from http://data-nctcoggis.opendata.arcgis.com/datasets/2015-land-use
north_texas_inputs <- list(
  name = 'North Texas', 
  path = '2015_North_Texas_Land_Use/2015_Land_Use.shp',
  feature = 'CATEGORY',
  sf = c('Single family'),
  mf= c('Multi-family'),
  or = c('Mixed use', 'Residential acreage'),
  separator = 'no_separator', 
  proj4string = lat_long
)


dallas_inputs <- list(
  name = 'Dallas',
  path = 'dallas_zoning/BaseZoning.shp',
  sf = c('A(A)', 'D(A)', 'TH', 'R'),
  mf = c('CH', 'MF', 'MU'),
  or = c('MH(A)'),
  separator = '-',
  feature = 'LONG_ZONE_',
  proj4string = texas_north_central
)

denton_inputs <- list(
  name = 'Denton',
  path ='denton_zoning/Current_zoning.shp',
  feature = 'ZONING',
  sf = c('DR-1', 'NR-1', 'NR-2', 'NR-3', 'NR-4', 'NR-5', 'NR-6', 'RD-5'),
  mf = c('DR-2', 'MF-1', 'NRMU', 'NRMU-12', 'RCR-1', 'RCR-2', 'RCC-D', 'RCC-N'),
  or = c('PD', 'RD-5X'),
  separator = 'no_separator', 
  proj4string = texas_north_central
)


if (getOption('run.main', default=TRUE)) {
  
  print('')
  

}
