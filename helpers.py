import numpy as np
import datetime as dt
import time
import pandas as pd
import geopandas as gpd
import fiona
import shapely.geometry
import matplotlib.pyplot as plt
from inputs import *

# Purpose of this file is to do non-spatial data processing. It begins with a basic helper function for color processing.
# Then, in this order, there are processing functions for (1) zip code boundaries for the entire US, (2) zoning data
# for Austin, Dallas, and Central North Texas, (3) Demand Data from Realtor and Zillow, (4) Parcel data for Austin and
# Dallas, (5) Permit data for Austin, Dallas, and Houston.

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


# Get zip boundaries data for the entire United States ----------------------------------------------------------------------- Zip Data
def get_zip_boundaries():
    zipdata = gpd.read_file(zip_boundaries_path)
    zipdata.index = zipdata['ZCTA5CE10']
    zipdata.index = [str(ind) for ind in zipdata.index]
    return zipdata

# Checks whether a string can be converted to a float
def will_it_float(text):
    try:
        text = float(text)
        return True
    except ValueError:
        return False

# Zoning data processing functions-------------------------------------------------------------------------------------------- Zoning Data

def process_zoning_shapefile(input, overlay = {'-MU':'Multifamily'}, broaden = True, parse_base_zones = True):
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

    if parse_base_zones and input.regulations_path is not None:

        # Get regulations data - we only read this in because we need the list of base zones.
        reg_data = pd.read_csv(input.regulations_path, index_col=0, encoding='Latin-1')

        # If the data includes non-base zones, only include base zones in this particular search
        if 'class' in reg_data.columns:
            reg_data = reg_data.loc[reg_data['class'] == 'base']

        reg_data.index = [str(ind) for ind in reg_data.index]

        # Process the zoning data - this function is well defined, I checked
        base_zones = reg_data.index.unique().tolist()

        def process_zone(text):
            for i in base_zones:
                if i in text:
                    return i
            return 'Unknown'

        # Apply
        raw_data['base_zone'] = raw_data[input.feature].apply(process_zone)

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

    # Helper functions
    # Convert to datetime
    def to_datetime(inputstring):
        astring = str(inputstring)
        year = int(astring.split('-')[0])
        month = int(astring.split('-')[1])
        return dt.date(year=year, month=month, day=1)

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

def process_dallas_parcel_data(quickly = False):
    """
    :param quickly: Boolean, if true, will read just the 2016 parcel data without processing it or joining it to other
    data. Otherwise will join data to 2013 parcel data.
    :return: Parcel data
    """

    if quickly:
        p2016 = gpd.read_file(dallas_parcel_data_path_2016)
        p2016 = p2016.rename(columns = {'sptbcode':'sptbcode_2016'})
        return p2016

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

def process_houston_parcel_data(feature_files = [harris_parcel_building_res_path_2018],
                                feature_columns_list = [houston_building_res_columns], county_level = False):
    """
    Merge subsetted houston parcel data with harris county data. By default uses a left join.
    :param feature_files: A list of paths of feature files to merge with the data. Defaults to
    the building_res file path. Can also be a string of one file (doesn't have to be a list). If you set this equal to
     None, it will just return the parcel shape data with no other data attached (although the parcel shape data does
     have a couple of useful factors).
    :param feature_columns_list: A list of column headers for each file - you can find these in the
     Harris_Parcel_Feature_Columns microsoft database.
     :param county_level: Boolean, default False. If true, will read in parcel data from all of Harris county, not just
     Houston.
    :return: GeoDataFrame with the parcel shapes and the merged data.
    """

    # Get parcel data
    if county_level:
        geodata = gpd.read_file(harris_parcel_path_2018)
    else:
        geodata = gpd.read_file(houston_parcel_path_2018)

    geodata = geodata.loc[geodata['geometry'].is_valid]
    geodata.reset_index(inplace = True)
    geodata = geodata.loc[geodata['HCAD_NUM'].apply(will_it_float)]
    geodata['HCAD_NUM'] = geodata['HCAD_NUM'].astype(float)

    # Make sure feature_files is a list not a string
    if isinstance(feature_files, str):
        feature_files = [feature_files]

    # Get data
    if feature_files is not None:
        for feature_path, feature_header in zip(feature_files, feature_columns_list):
            data = pd.read_csv(feature_path, sep = '\t', header = None, encoding='Latin-1')
            data.columns = feature_header
            data = data.rename(columns = {'ACCOUNT':'HCAD_NUM'})
            data = data.drop_duplicates(subset=['HCAD_NUM'], keep='first')
            print('Starting to merge, geodata shape is {}, data shape is {}'.format(geodata.shape, data.shape))
            geodata = geodata.merge(data, how = "left", on = "HCAD_NUM")

    # Return
    return geodata


# Permit data -------------------------------------------------------------------------------------------------------------- permit data

def process_austin_permit_data(searchfor, permittypedesc = None, workclass = None,
                               earliest = None, latest = None):
    """
    :param searchfor: List of strings to search for in the 'PermitClass' and 'Description' columns. Often it's worth
    explicitly passing in permit classes, i.e. searchfor = ['101 single family houses'].
    :param permittypedesc: The permittypedesc to match. Ex: "Building Permit."
    :param workclass: Workclass to match. Ex: "New"
    :param earliest: Earliest date to consider (inclusive). Data runs from 1971-2018.
    :param latest: Latest date to consider (inclusive). Data runs from 1971-2018.
    :return: GeoDataFrame of subsetted permit data
    """

    time0 = time.time()
    print('Reading Austin permit data')
    permit_data = pd.read_csv(austin_permit_path)
    print('Finished reading Austin permit data, took {}. Now subsetting.'.format(time.time() - time0))

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

    print('Finished processing Austin permit data, took {}.'.format(time.time() - time0))
    return construction

# Pretty self explanatory, except ONLY use this if you don't want to use the corrected data. Dallas permit data
# has some incorrect lat and longs near the center of the city, which have been corrected and stored at
# dpm_save_path (listed in inputs) using the zipcode_scrape.py functions.
def process_dallas_permit_data(permit_types, earliest = None, latest = None):
    """
    :param permit_types: List of permit types to filter for. Will only consider rows where the permit type is one of these
    permit types.
    :param earliest: Earliest date to consider. Data runs from 2011-2016. Defaults to None.
    :param latest: Latest date to consider. Data runs from 2011-2016. Defaults to None.
    :return: GeoDataFrame of subsetted permit data
    """

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


def process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE'], 
                                searchin = ['PROJ_DESC'], 
                                kind = 'structural',
                                earliest = None, latest = None):
    """
    Process houston permit data. Note that houston permit data does not specify whether housing is new or not, so I have
    to parse the descriptions to figure out whether housing is new sf/new mf housing. Thankfully the descriptions
    are rather formulaic so this is not too hard to do.
    :param searchfor: A list of keywords to search for, i.e. ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE']).
    To subset to only include new construction, each string should start with 'NEW '.
    return
    :param searchin: The columns to search in. Defualts to ['PROJ_DESC'], project description. Note that the function
    will return rows where ANY of the columns specified by 'searchin' contain ANY of the keywords in 'searchfor'.
    :param kind: The kind of permit data to read in. Defaults to 'structural', can either by 'structural' 
    or 'demolition'
    :param earliest: Earliest date to consider.
    :param latest: Latest date to consider.
    :return: GeoDataFrame of subsetted permit data.
    """
    # Load data
    if kind == 'structural':
        path = houston_structural_permits_path
    elif kind == 'demolition':
        path = houston_demolition_permits_path
    else:
        print("""Error: in process_houston_permit_data call, "kind" arg must either 
            equal "demolition" or "structural", not "{}""".format(kind))
        return None
    permit_data = gpd.read_file(path)
    permit_data = permit_data.rename(columns = {'SITUS_ZIP_':'Zipcode'}) # Already has 100% complete zipcode information.
    permit_data = permit_data.loc[~permit_data['geometry'].apply(lambda x: x is None)]

    # Ignore repeated plans
    permit_data = permit_data.loc[~permit_data['PROJ_DESC'].str.contains('REPEAT PLAN')]

    # Search for kwargs in the specified columns, almost always will be in project description
    flags = pd.Series(False, index = permit_data.index)
    for col in searchin:
        flags = flags | permit_data[col].str.contains('|'.join(searchfor))
    subsetted_data = permit_data.loc[flags]

    # Process dates. Note that data runs from 1988 to 2018.
    subsetted_data.loc[:, 'Year'] = subsetted_data['APPLN_DATE'].apply(lambda x: int(x[0:4]))
    if earliest is not None:
        subsetted_data = subsetted_data.loc[subsetted_data['Year'] >= earliest]
    if latest is not None:
        subsetted_data = subsetted_data.loc[subsetted_data['Year'] <= latest]

    # Process whether permit has been approved yet (using scraped data)
    approvals = pd.read_csv('shared_data/houston_permit_statuses.csv', index_col = 0)
    approvals.index = [str(id) for id in approvals.index]

    # Merge
    subsetted_data['Approval'] = subsetted_data['PROJECT_NO'].map(approvals['Approval'])
        
    # Return
    return subsetted_data

