import numpy as np
import datetime as dt
import time
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from inputs import *

# To do: Need to simplify very complex shapefiles better than geopandas does it.

# Very basic helper functions ---------------------------------------------------------------------------------------

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

# Zoning data processing functions------------------------------------------------------------------------------------

def process_zoning_shapefile(input, broaden = True):
    """"
    Processes shapefile.
    :param input: An input, of class 'zoning_inputs.'
    :param broaden: Boolean. Default true. If true, will decode zoning data into a broader classification (i.e. sf, mf)
    as specified by the zoning input class - this processed zoning data will go in a column labelled "broad_zone."
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

    def get_zone(text):
        split_text = text.split(input.separator)[0]
        for key in input.base_zones:
            if split_text in input.base_zones[key]:
                return(key)
        return([key for key in input.base_zones][-1])

    raw_data['broad_zone'] = raw_data[input.feature].apply(get_zone)

    # Manually set CRS if necessary
    if input.crs is not None:
        raw_data.crs = input.crs

    print('Finished processing zones and crs, took {}'.format(time.time() - time0))


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

# Demand data processing functions, sometimes with limited graphing ability -------------------------------------------
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

# Get zip boundaries data
def get_zip_boundaries():
    zipdata = gpd.read_file(zip_boundaries_path)
    zipdata.index = zipdata['ZCTA5CE10']
    zipdata.index = [str(ind) for ind in zipdata.index]
    return zipdata


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

    # This doesn't work
    #print('Finished merging files, took {}. Now saving file.'.format(time.time() - time0))
    #parcel_data.to_file(dallas_processed_parcel_data_path)
    #print('Finished, took {}'.format(time.time() - time0))

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



if __name__ == '__main__':

    pass