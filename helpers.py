import numpy as np
import datetime as dt
import time
import pandas as pd
import geopandas as gpd
import fiona
import shapely.geometry
import matplotlib.pyplot as plt
from inputs import *

# To do: Need to simplify very complex shapefiles better than geopandas does it.

# Very basic helper functions ---------------------------------------------------------------------------------------------- Helper functions

# Helper color function - adapted from https://stackoverflow.com/questions/35516318/plot-colored-polygons-with-geodataframe-in-folium
def convert_to_hex(rgba_color):
    red = str(hex(int(rgba_color[0]*255)))[2:].capitalize()
    green = str(hex(int(rgba_color[1]*255)))[2:].capitalize()
    blue = str(hex(int(rgba_color[2]*255)))[2:].capitalize()

    if blue=='0':
        blue = '00'
    if red=='0':
        red = '00'
    if green=='0':
        green='00'

    return '#'+ red + green + blue

# Convert to datetime
def to_datetime(inputstring):
    astring = str(inputstring)
    year = int(astring.split('-')[0])
    month = int(astring.split('-')[1])
    return dt.date(year = year, month = month, day = 1)

def quiet_to_datetime(inputstring):
    try:
        result = to_datetime(inputstring)
        return result
    except:
        return inputstring

# Check which row/column entries are datetimes,
# and put the rest in the multiindex
def in_multiindex(astring):
    try:
        to_datetime(astring)
        return False
    except:
        return True


# Get zip boundaries data
def get_zip_boundaries():
    zipdata = gpd.read_file(zip_boundaries_path)
    zipdata.index = zipdata['ZCTA5CE10']
    zipdata.index = [str(ind) for ind in zipdata.index]
    return zipdata

# Zoning data processing functions-------------------------------------------------------------------------------------------- Zoning Data

def process_zoning_shapefile(input, broaden = True, overlay = {'-MU':'Multifamily'}):
    """"
    Processes shapefile.
    :param input: An input, of class 'zoning_inputs.'
    :param broaden: Boolean. Default true. If true, will decode zoning data into a broader classification (i.e. sf, mf)
    as specified by the zoning input class - this processed zoning data will go in a column labelled "broad_zone."
    :param overlay: A dictionary which maps overlay strings to broadened outputs. This is mostly useful for dealing with
    mixed use overlays in austin (there are no mixed-use base zones, so passing in {'-mu':'Multifamily'} will help the
    function recognize that mixed-use zones are multifamily zones.
    """

    # Begin timing
    time0 = time.time()

    # Make sure input is correct
    if isinstance(input, zoning_inputs) == False:
        warning('Error, Input must be of class zoning_inputs')
        return None

    # Read data and process zone codes
    print('Reading file')
    raw_data = gpd.read_file(input.path)
    print('Finished reading file, took {}'.format(time.time() - time0))

    if broaden:
        def get_zone(text):
            split_text = text.split(input.separator)[0]
            for key in input.base_zones:
                if split_text in input.base_zones[key]:
                    return(key)
            return([key for key in input.base_zones][-1])

        raw_data['broad_zone'] = raw_data[input.feature].apply(get_zone)

        # Now put account for overlays
        if overlay is not None:
            for key in [key for key in overlay]:
                raw_data.loc[raw_data[input.feature].str.contains(key), 'broad_zone'] = overlay[key]

    # Manually set CRS if necessary
    if input.crs is not None:
        raw_data.crs = input.crs

    print('Finished processing zones, mixed use overlay, and crs, took {}'.format(time.time() - time0))


    return raw_data

def get_regulation_data(zones, regulations_path, regulation_feature, fill = None):
    """
    :param zones: a pandas series of zoning.
    :param regulations_path: The path for the data which specifies regulations by base zone or zone code.
    :param regulation_feature: The feature of interest to find by zone.
    :param fill: How to fill na values inside the specific regulation feature, i.e. 0. Defaults to None, will leave empty
    values as empty.
    :return: A series of the desired regulation feature, i.e. minimum lot size.
    """

    # Get regulations data
    reg_data = pd.read_csv(regulations_path, index_col = 0, encoding = 'Latin-1')

    # If the data includes non-base zones, only include base zones in this particular search
    if 'class' in reg_data.columns:
        reg_data = reg_data.loc[reg_data['class'] == 'base']

    reg_data = reg_data[regulation_feature]
    reg_data.index = [str(ind) for ind in reg_data.index]
    # Add missing object just in case


    # Fill missing values
    if fill is None:
        pass
    else:
        reg_data = reg_data.fillna(fill)
        reg_data['Unknown'] = fill

    # Process the zoning data - this function is well defined, I checked
    base_zones = reg_data.index.unique().tolist()
    def process_zone(text):
        for i in base_zones:
            if i in text:
                return i
        return 'Unknown'

    # Process and then you're done
    pzones = zones.map(process_zone)
    result = pzones.map(reg_data)
    return result

# Demand data processing functions, sometimes with limited graphing ability -------------------------------------------------- Demand Data
def process_demand_data(input, graph = False, style = 'Line', date = dt.date(year = 2018, month = 4, day = 1)):
    """
    Initially process and graph housing demand datasets. If it graphs, it will not be a choropleth.
    :param input: A Demand_Input class.
    :return: data (a pd dataframe of just the time series and the index),
    metadata (a pd dataframe of metadata about the index)
    """

    # Check that it's the right class
    if isinstance(input, Demand_Input) == False:
        warning('Error, Input must be of class Demand_Input')
        return None

    path = input.path
    geo_filter = input.geo_filter
    geo_filter_values = input.geo_filter_values
    index_col = input.index_col

    # Read data
    raw_data = pd.read_csv(path, index_col = index_col)

    # Only consider specific geographic subset
    raw_data = raw_data.loc[raw_data[geo_filter].isin(geo_filter_values), :]

    # Separate metadata and data
    if input.source == 'Zillow':
        data_cols = raw_data.columns[[not in_multiindex(astring) for astring in raw_data.columns]]
        metadata_cols = raw_data.columns[[in_multiindex(astring) for astring in raw_data.columns]]
        metadata = raw_data[metadata_cols]
        data = raw_data[data_cols].transpose()
        data.index = [quiet_to_datetime(astring) for astring in data.index]
    elif input.source == 'Realtor':
        data = raw_data[[input.feature]]
        metadata = raw_data[[geo_filter]]
        data.index = [str(int(ind)) for ind in data.index]
        metadata.index = [str(int(ind)) for ind in metadata.index]

    # Now graph
    if graph:

        xtick_rotation = 90
        xtick_size = 8

        if input.source == 'Zillow' and style == 'Line':
            for location in geo_filter_values:
                print('Starting to graph for {}'.format(location))
                ids = [region for region in data.columns if metadata.loc[region, geo_filter] == location]
                filtered_data = data[ids]
                filtered_data.plot(title = location + ' ,' + input.name)
                plt.show()

        elif input.source == 'Zillow':
            for location in geo_filter_values:
                print('Starting to graph for {}'.format(location))
                for region in data.columns:
                    print(metadata.loc[region, geo_filter])
                ids = [region for region in data.columns if metadata.loc[region, geo_filter] == location]
                filtered_data = data.loc[date, ids].sort_values()
                filtered_data.plot(kind = 'bar', color = 'Green', title = location + ', ' + input.name, legend = False)
                plt.xticks(rotation=xtick_rotation, size = xtick_size)
                plt.show()

        elif input.source == 'Realtor':
            for location in input.geo_filter_classes:
                print('Starting to graph for {}'.format(location))
                ids = [region for region in data.index if metadata.loc[region, geo_filter] == location]
                filtered_data = data.loc[ids].sort_values(input.name)
                filtered_data.plot(kind = 'bar', title = location)
                plt.xticks(rotation=xtick_rotation, size = xtick_size)
                plt.show()

    return data, metadata

# Add demand noncategorical data to zip geodata
def add_demand_data(zip_geodata, demand_input, city, feature_name = None):
    """
    :param zip_geodata: Zip geodata, should already have all the geometries for all of the zip codes in the city.
    :param demand_input: A class demand_input.
    :param cities: Iterable of cities, i.e. ['Austin, TX'] or ['Austin, TX', 'Dallas, TX']
    :param feature_name: Updated feature name in case of conflicting features from different sources (i.e. 'Avg Listing Price' is in both the sfhomes and the cth homes dataset)
    :return: zip_geodata updated with the new value
    """

    if feature_name is None:
        feature_name = demand_input.feature

    data, metadata = process_demand_data(demand_input, graph = False)
    data = data.loc[[index for index in data.index if metadata.loc[index, demand_input.geo_filter] == city]]
    zip_geodata[feature_name] = data[demand_input.feature]
    return zip_geodata


# Parcel data -------------------------------------------------------------------------------------------------------------- parcel data

# Join old and new dallas parcel data
def join_dallas_parcel_data():

    time0 = time.time()
    print('Reading 2013 file')
    p2013 = gpd.read_file(dallas_parcel_data_path_2013)
    p2013.drop_duplicates('gis_acct', inplace = True) # Only drops about 700/300K
    p2013.columns = [str(col) + '_2013' for col in p2013.columns]
    p2013 = p2013.rename(columns = {'gis_acct_2013': 'gis_acct'})

    print('Finished reading 2013 file, took {}. Now reading 2016 file.'.format(time.time() - time0))
    p2016 = gpd.read_file(dallas_parcel_data_path_2016)
    p2016.drop_duplicates('gis_acct', inplace = True) # Only drops about 750/400K
    p2016.columns = [str(col) + '_2016' for col in p2016.columns]
    p2016 = p2016.rename(columns = {'gis_acct_2016': 'gis_acct'})


    print('Finished reading 2016 file, took {}. Now merging files.'.format(time.time() - time0))
    parcel_data = p2016.merge(p2013, on = 'gis_acct')
    parcel_data = gpd.GeoDataFrame(data = parcel_data[[col for col in parcel_data.columns if 'geometry' not in col]],
                                   geometry = parcel_data['geometry_2016'])
    print(parcel_data.columns, parcel_data.index)
    return parcel_data

# Permit data -------------------------------------------------------------------------------------------------------------- permit data

def process_dallas_parcel_data():

    time0 = time.time()

    print('Loading data')
    parcel_data = join_dallas_parcel_data()
    print('Finished loading data, took {} seconds. Now normalizing data conditionally on land use.'.format(time.time() - time0))

    parcel_data['area'] = 1000*parcel_data['geometry'].area

    parcel_data.loc[:, 'value_per_area'] = parcel_data['tot_val_2013'].divide(parcel_data['area'])
    land_use_means = {}
    land_use_stds = {}
    for code in parcel_data['sptbcode_2013'].unique():
        subset = parcel_data.loc[parcel_data['sptbcode_2013'] == code]
        land_use_means[code] = subset['value_per_area'].mean()
        land_use_stds[code] = subset['value_per_area'].std()

    # Normalize - there should be no key errors because we just called parcel_data['land_use'].unique().
    vpa_mean = parcel_data['sptbcode_2013'].map(land_use_means)
    vpa_std = parcel_data['sptbcode_2013'].map(land_use_stds)
    norm_pva = (parcel_data['value_per_area'] - vpa_mean).divide(vpa_std)
    parcel_data['normalized_value_per_area'] = norm_pva

    return parcel_data

def process_austin_permit_data(searchfor, permittypedesc = None, workclass = None,
                               earliest = None, latest = None):
    """
    :param searchfor:
    :param permittypedesc: The permittypedesc to match. Ex: "Building Permit."
    :param workclass: Workclass to match. Ex: "New"
    :param earliest:
    :param latest:
    :return:
    """

    time0 = time.time()
    print('Reading permit data')
    permit_data = pd.read_csv(austin_permit_path)
    print('Finished reading permit data, took {}. Now subsetting.'.format(time.time() - time0))

    # Initial subset
    construction = permit_data
    if permittypedesc is not None:
        construction = construction.loc[(permit_data['PermitTypeDesc'] == permittypedesc)]
    if workclass is not None:
        construction = construction.loc[(permit_data['WorkClass'] == workclass)]

    # Get right timeframe
    if earliest is not None:
        construction = construction.loc[construction['CalendarYearIssued'] >= earliest, :]

    if latest is not None:
        construction = construction.loc[construction['CalendarYearIssued'] <= latest, :]

    def to_lowercase(object):
        return str(object).lower()

    # Now subset to just include entries where one of the 'searchfor' strings is in either the description or permit
    #  class column. Can use this to search for either multifamily or single family, or anything specific really.
    # Best to search based on permit class, see below.
    searchfor = [str(x).lower() for x in searchfor]
    construction.loc[:, 'PermitClass'] = construction.loc[:, 'PermitClass'].apply(to_lowercase)
    construction.loc[:, 'Description'] = construction.loc[:, 'Description'].apply(to_lowercase)
    construction = construction.loc[(construction['PermitClass'].str.contains('|'.join(searchfor))) |
                                    (construction['Description'].str.contains('|'.join(searchfor))), :]

    def return_point(series):
        # series[0] is long, series[1] is lat
        return shapely.geometry.point.Point(series[0], series[1])

    # Now process geometry
    points = construction[['Longitude', 'Latitude']].apply(return_point, axis = 1)
    construction = gpd.GeoDataFrame(data = construction, geometry = points)
    construction.crs = {'init': 'epsg:4326'}

    print('Finished processing permit data, took {}.'.format(time.time() - time0))
    return construction

# For Austin, calendar year data
# np.sort(construction['CalendarYearIssued'].unique())
# array([1971, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986,
       #1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997,
       #1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
       #2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018],
      #dtype=int64)
# For Austin, list of construction['PermitClass'].unique() is below (took the liberty of sorting)
# array(['c - medical gas', 'c- 101 single family houses',
       #'c- 103 two family bldgs', 'c- 104 three & four family bldgs',
       #'c- 105 five or more family bldgs', 'c- 106 mixed use',
       #'c- 213 hotels, motels, & tourist cabins',
       #'c- 214 other nonhousekeeping shelter',
       #'c- 318 amusement, social & rec bldgs',
       #'c- 319 churches and othr religious bldgs',
       #'c- 320 industrial bldgs', 'c- 321 pkg garage bldg & open deck',
       #'c- 322 service station & repair garage',
       #'c- 323 hospital & institutional bldgs',
       #'c- 324 office, bank & professional bldgs',
       #'c- 325 public works & utilities bldgs',
       #'c- 326 schools & other educational bldgs',
       #'c- 327 stores & customer services',
       #'c- 328 commercial other nonresident bldg',
       #'c- 329 com structures other than bldg',
       #'c- 437 addn, alter, convn-nonres',
       #'c- 647 demolition 3 and 4 family bldgs',
       #'c- 648 demolition 5 or more family bldgs',
       #'c- 649 demolition all other bldgs com',
       #'c-1000 commercial remodel', 'c-1001 commercial finish out',
       #'c-1002 commercial remodel & finish out',
       #'c-2000 relocation commercial', 'c-2002 commercial boat dock',
       #'c-2003 commercial retaining wall', 'com. curb, gutter',
       #'com. driveway', 'com. driveway & sidewalk',
       #'com. driveway, curb, gutter',
       #'com. driveway,sidewalk,curb,gutter', 'com. sidewalk',
       #'com. sidewalk, curb, gutter', 'commercial', 'nan',
       #'r- 101 single family houses', 'r- 102 secondary apartment',
       #'r- 103 two family bldgs',
       #'r- 328 resident other nonresident bldg',
       #'r- 329 res structures other than bldg',
       #'r- 330  accessory use to primary',
       #'r- 434 addition & alterations', 'r- 435 renovations/remodel',
       #'r- 436 addn to increase housing units',
       #'r- 437 residential boat dock',
       #'r- 438 residential garage/carport addn',
       #'r- 438 residential retaining wall',
       #'r- 645 demolition one family homes',
       #'r- 646 demolition two family bldgs',
       #'r- 649 demolition all other bldgs res',
       #'r-2001 relocation residential', 'replace of curb & gutter',
       #'res. curb, gutter', 'res. driveway', 'res. driveway & sidewalk',
       #'res. driveway, curb, gutter',
       #'res. driveway, sidewalk, curb, gutter', 'res. sidewalk',
       #'res. sidewalk, curb, gutter', 'residential', 'sign permit'],
      #dtype=object)

# Columns: Index(['PermitType', 'PermitTypeDesc', 'PermitNum', 'PermitClassMapped',
       #'PermitClass', 'WorkClass', 'Condominium', 'ProjectName', 'Description',
       #'TCAD_ID', 'PropertyLegalDescription', 'AppliedDate', 'IssuedDate',
       #'DayIssued', 'CalendarYearIssued', 'FiscalYearIssued',
       #'IssuedInLast30Days', 'IssuanceMethod', 'StatusCurrent', 'StatusDate',
       #'ExpiresDate', 'CompletedDate', 'TotalExistingBldgSQFT',
       #'RemodelRepairSQFT', 'TotalNewAddSQFT', 'TotalValuationRemodel',
       #'TotalJobValuation', 'NumberOfFloors', 'HousingUnits',
       #'BuildingValuation', 'BuildingValuationRemodel', 'ElectricalValuation',
       #'ElectricalValuationRemodel', 'MechanicalValuation',
       #'MechanicalValuationRemodel', 'PlumbingValuation',
       #'PlumbingValuationRemodel', 'MedGasValuation', 'MedGasValuationRemodel',
       #'OriginalAddress1', 'OriginalCity', 'OriginalState', 'OriginalZip',
       #'CouncilDistrict', 'Jurisdiction', 'Link', 'ProjectID',
       #'MasterPermitNum', 'Latitude', 'Longitude', 'Location',
       #'ContractorTrade', 'ContractorCompanyName', 'ContractorFullName',
       #'ContractorPhone', 'ContractorAddress1', 'ContractorAddress2',
       #'ContractorCity', 'ContractorZip', 'ApplicantFullName',
       #'ApplicantOrganization', 'ApplicantPhone', 'ApplicantAddress1',
       #'ApplicantAddress2', 'ApplicantCity', 'ApplicantZip', 'geometry'],
      #dtype='object')


# Pretty self explanatory, except ONLY use this if you don't want to use the corrected data. Dallas permit data
# has some incorrect lat and longs near the center of the city, which have been corrected and stored at
# dpm_save_path (listed in inputs) using the zipcode_scrape.py functions.
def process_dallas_permit_data(permit_types, earliest = None, latest = None):

    permit_data = pd.read_csv(dallas_permit_path)

    # Get date
    def get_year(text):
        date = text.split(' ')[0]
        year = date.split('/')[-1]
        return int(year)

    permit_data['Year'] = permit_data['Issued Date'].apply(get_year)

    # subset by year and permit type, then return
    construction = permit_data.loc[permit_data['Permit Type'].isin(permit_types)]
    if earliest is not None:
        construction = construction.loc[construction['Year'] >= earliest]
    if latest is not None:
        construction = construction.loc[construction['Year'] <= latest]

    # Get points
    def return_point(location):
        coords = str(location).split('\n')[-1]
        lat = float(coords[1:].split(', ')[0])
        long = float(coords.split(', ')[-1][:-1])
        return shapely.geometry.point.Point(long, lat)

    # Geocoding fails for 5000/220,000 locations, which have addresses but not lat/long. Could improve this using google API
    # or using a dictionary of addresses:lat/long (because some addresses are repeated) but a bit pressed for time right now.
    # Fails for about 15% of construction sites.
    def is_geocoded(location):
        try:
            return_point(location)
            return True
        except:
            return False

    construction = construction.loc[construction['GeoLocation'].apply(is_geocoded)]
    points = construction['GeoLocation'].apply(return_point)
    construction = gpd.GeoDataFrame(data = construction, geometry = points)

    return construction