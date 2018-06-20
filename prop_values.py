import os
import sys
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import shapefile_work_around as swa
import geopandas as gpd
import folium
import branca

os.chdir("C:/Users/amspe/Documents/R/MI2018/TXHousing")

# Global paths
zip_boundaries_path = "data\cb_2017_us_zcta510_500k\cb_2017_us_zcta510_500k.shp"
nbhd_boundaries_path = "data\Zillow Data\ZillowNeighborhoods-TX\ZillowNeighborhoods-TX.shp"

# Helper functions -------------------------------------------------------------------------------------------

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

# Classes --------------------------------------------------------------------------------------------------------
class City:
    def __init__(self, name, long, lat, zoom, xlims = [None, None], ylims = [None, None]):
        self.name = name
        self.long = long
        self.lat = lat
        self.zoom = zoom
        self.xlims = xlims
        self.ylims = ylims

# Cities (needed for defining Demand_Input)
austin = City('Austin', long=30.267, lat=-97.743, zoom=11, xlims = [29.25, 31.75], ylims = [96.5, 99])
dallas = City('Dallas', long = 32.7767, lat=-96.7970, zoom = 11, xlims = [-97.3, -96.1], ylims = [32.25, 33.25])
houston = City('Houston', long = 29.7604, lat = -95.3698, zoom = 11, xlims = [28, 31], ylims = [-93, -97])

# Input format
class Demand_Input:

    def __init__(self, path, source,
                 feature = 'Value',
                 geo_filter_classes=[austin, dallas, houston],
                 geography = 'zip',
                 **kwargs):
        """
        :param path: Path of the dataset, assumed csv
        :param source: Source of the data. Can either be "Zillow" or "Realtor"
        :param feature: Name of the feature. In the case of Realtor data, this must be the column of the data.
        :param geo_filter_classes: A list of City classes which are used to subset the data.
        :param geography: The level of detail of the dataset, i.e. "zip" or "neighborhood."
        :param geo_filter: A kwarg. This is the geography level (i.e. city, state, county) by which to filter data,
        and it should be a column of the dataset. Default depends on the source.
        :param index_col: A kwarg. The column to use as the index. Default depends on the source.
        :param name: an alternate name for graphical display of the feature. Defaults to the feature string.
        """


        self.path = path
        self.source = source
        self.feature = feature
        self.geo_filter_classes = geo_filter_classes
        self.geography = geography

        # Depending on the source, use different defaults
        try:
            self.geo_filter == kwargs['geo_filter']
        except:
            if source == 'Zillow':
                self.geo_filter = 'City'
            elif source == 'Realtor':
                self.geo_filter = 'ZipName'

        try:
            self.index_col = kwargs['index_col']
        except:
            if source == 'Zillow':
                self.index_col = 'RegionName'
            elif source == 'Realtor':
                self.index_col = 'ZipCode'


        # Get filter values and name from kwargs
        try:
            self.geo_filter_values = kwargs['geo_filter_values']
        except:
            self.geo_filter_values = [city.name for city in self.geo_filter_classes]

        try:
            self.name = kwargs['name']
        except:
            self.name = feature

# Inputs ----------------------------------------------------------------------------------------------------------

# Zillow inputs
sfhomes_nbhd = Demand_Input(path = 'data/Zillow Data/sfhomes_neighborhood.csv',
                        source = 'Zillow',
                        feature = 'Single Family Homes Median Value',
                        geography = 'Neighborhood')

inventoryraw_zip = Demand_Input(path = "data/Zillow Data/inventoryraw_zip.csv",
                            source = 'Zillow',
                            feature = 'Raw Inventory',
                            geography = 'Zip')

medpricecuts_zip = Demand_Input(path = "data/Zillow Data/medpricecut_zip.csv",
                            source = 'Zillow')

allhomeprices_zip = Demand_Input(path = 'data/Zillow Data/allhomeprice_zip.csv',
                             source = 'Zillow')

# One note is that to avoid a rather annoying index error I had to change the name of "Downtown" region in Austin to
# "downtown" (without the capital D).
allhomeprices_nbhd = Demand_Input(path = "data/Zillow Data/allhomeprice_neighborhood.csv",
                                  source = 'Zillow',
                                  geography = 'Neighborhood')

# Realtor inputs
realtor_med_dom_cbsa = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Median DOM (vs CBSA)',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_med_dom_yoy = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Median DOM Y/Y',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_med_vpp_cbsa = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Views Per Property  (vs CBSA)',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])


realtor_med_vpp_yoy = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Views Per Property Y/Y',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_hotness_cbsa = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Hotness Rank Within CBSA',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_hotness_yoy = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Hotness Rank Y/Y',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

# Core realtor data
realtor_avg_price = Demand_Input(path = 'data/RDC_InventoryCoreMetrics_Zip.csv',
                                 source = 'Realtor',
                                 feature = 'Avg Listing Price',
                                 geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_avg_sf_price = Demand_Input(path = "data/RDC_InventoryCoreMetrics_Zip_sfh.csv",
                                 source = 'Realtor',
                                 feature = 'Avg Listing Price',
                                 geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'],
                                 name='Avg Listing Price, SF Homes')


realtor_avg_cth_price = Demand_Input(path = "data/RDC_InventoryCoreMetrics_Zip_cth.csv",
                                 source = 'Realtor',
                                 feature = 'Avg Listing Price',
                                 geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'],
                                 name='Avg Listing Price, Condo / Townhomes')



# Graphing and processing functions -------------------------------------------------------------------------------

def process_and_graph(input, graph = False, style = 'Line', date = dt.date(year = 2018, month = 4, day = 1)):
    """
    Initially process and graph housing datasets.
    :param input: A Zillow_Inputs class.
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

def choropleth_processed_data(input, date = dt.date(year = 2018, month = 4, day = 1), plot_folium = False, colormap = 'Blues'):
    """
    :param input: Input class
    :param date: Date of data to graph
    :param plot_folium: If true, uses folium (this is a good idea)
    :param geography: Default 'zip'. Otherwise, will try to plot data at the neighborhood level.
    :param kwargs: Other graphing kwargs
    :return:
    """

    # Get data and metadata
    data, metadata = process_and_graph(input, graph = False)

    # Get shapefile, reading this in takes 10 sec or so
    if input.geography == 'zip':
        boundaries = gpd.read_file(zip_boundaries_path)
        boundaries.index = boundaries['ZCTA5CE10']
    else:
        boundaries = gpd.read_file(nbhd_boundaries_path)
        boundaries.index = boundaries['Name']

    # Merge
    if input.source == 'Zillow':
        data.columns = [str(code) for code in data.columns]
        metadata.index = [str(index) for index in metadata.index]
        boundaries[input.name] = data.loc[date]
    elif input.source == 'Realtor':
        data.index = [str(int(code)) for code in data.index]
        metadata.index = [str(int(code)) for code in metadata.index]
        boundaries[input.name] = data[input.feature]

    boundaries = boundaries.loc[boundaries[input.name].notnull()]

    # Plot
    if plot_folium:
        boundaries['index'] = boundaries.index
        for location in input.geo_filter_classes:
            basemap = folium.Map([location.long, location.lat], zoom_start=location.zoom)
            basemap.choropleth(geo_data=boundaries,
                               data=boundaries,
                               columns=['index', input.name],
                               key_on='feature.properties.index',
                               fill_color=colormap,
                               legend_name=input.name,
                               fill_opacity=0.8,
                               line_opacity=0)
            basemap.save('Figures/Testing/' + location.name + 'map.html')

    else:
        # Plot by location
        for location in input.geo_filter_values:
            print('Starting to graph for {}'.format(location))

            if input.source == 'Zillow':
                ids = [region for region in data.columns if metadata.loc[region, input.geo_filter] == location]
                filtered_data = boundaries.loc[ids]

            elif input.source == 'Realtor':
                ids = [region for region in data.index if metadata.loc[region, input.geo_filter] == location]
                filtered_data = boundaries.loc[ids].sort_values(input.name)

            filtered_data.plot(column = input.name, legend = True, cmap = colormap)
            plt.title(location + ', ' + input.name)
            plt.show()

#process_and_graph(sfhomes_nbhd, graph = True, style = 'bar')
choropleth_processed_data(realtor_avg_cth_price, plot_folium=True)