"""Processes property data, supplied by Realtor and Zillow"""

import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

# Helper functions ------------------------------------------------------

# Convert to datetime
def to_datetime(inputstring):
    astring = str(inputstring)
    year = int(astring.split('-')[0])
    month = int(astring.split('-')[1])
    return dt.date(year=year, month=month, day=1)

# Convert to datetime if possible
def quiet_to_datetime(inputstring):
    try:
        result = to_datetime(inputstring)
        return result
    except:
        return inputstring

# Check which row/column entries are datetimes,
# and put the rest in the multiindex
def in_multiindex(inputstring):
    try:
        to_datetime(inputstring)
        return False
    except:
        return True

# Class for property data ---------------------------------------------------------------------------------------------
class Property_Input():
    """
    :param path: Path of the dataset of property data. Data should be in a csv.
    :param source: Source of the data; can either be "Zillow" or "Realtor"
    :param feature: Name of the feature. In the case of Realtor data, this must be the column of the data.
    :param geo_filter: A kwarg. This is the geography level (i.e. city, state, county) by which to filter data,
        and it should be a column of the dataset. Default depends on the source.
    :param geo_filter_values: A list of values to subset the data to (by geo_filter). Defaults to some flavor of
        ['Austin', 'Houston', 'Dallas'] - but the specific strings are set to the source. If geo_filter_values == 'all',
        then the data will not be filtered.
    :param index_col: A kwarg. The column to use as the index. Default depends on the source.
    """

    def __init__(self, path, source, geo_filter = None, geo_filter_values = None, index_col = None):

        self.path = path
        if source not in ['Zillow', 'Realtor']:
            TypeError('source must equal either "Zillow" or "Realtor"')
        self.source = source

        # Depending on the source, use different defaults for geo_filter, geo_filter_values, index_col, and name
        if geo_filter is None:
            if source == 'Zillow':
                self.geo_filter = 'City'
            elif source == 'Realtor':
                self.geo_filter = 'ZipName'

        if geo_filter_values is None:
            if self.source == 'Zillow':
                self.geo_filter_values = ['Austin', 'Dallas', 'Houston']
            elif self.source == 'Realtor':
                self.geo_filter_values = ['Austin, TX', 'Dallas, TX', 'Houston, TX']
        else:
            self.geo_filter_values = geo_filter_values

        if index_col is None:
            if source == 'Zillow':
                self.index_col = 'RegionName'
            elif source == 'Realtor':
                self.index_col = 'ZipCode'


    def process_property_data(self, features = None, geo_filter_values = None):
        """
        :param features: For realtor data, subset to only include these features. Defaults to None.
        :param geo_filter_values: For convenience, you have the opportunity to override the self.geo_filter_values value
            with a new value, e.g. 'all' if you don't want to subset the data.
        :return: Two pandas dataframes, data and metadata (in that order). For Zillow inputs, data will have geographies
            in the index and a time series on the column, and metadata will have a variety of information about zip codes.
            For Realtor data, the data will have geographies in the index and a variety of features in the columns.

        Examples::

            realtor_hotness_data = Property_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv", source = 'Realtor')
            data, metadata = realtor_hotness_data.process_property_data(features = ['Median DOM (vs CBSA)', 'Views Per Property  (vs CBSA)']

        """
        if isinstance(features, str):
            features = [features]

        # Read data
        raw_data = pd.read_csv(self.path, index_col = self.index_col)

        # Only consider specific geographic subset
        if self.geo_filter_values != 'all':
            raw_data = raw_data.loc[raw_data[self.geo_filter].isin(self.geo_filter_values), :]

        # Separate metadata and data
        if self.source == 'Zillow':
            data_cols = raw_data.columns[[not in_multiindex(astring) for astring in raw_data.columns]]
            metadata_cols = raw_data.columns[[in_multiindex(astring) for astring in raw_data.columns]]
            metadata = raw_data[metadata_cols]
            data = raw_data[data_cols].transpose()
            data.index = [quiet_to_datetime(astring) for astring in data.index]
        elif self.source == 'Realtor':

            # Get rid of NaNs in index
            raw_data = raw_data.loc[raw_data.index.notnull()]

            # Subset if desired
            if features is not None:
                features.append(self.geo_filter)
                data = raw_data[features]
            else:
                data = raw_data
            metadata = raw_data[[self.geo_filter]]
            data.index = [str(int(ind)) for ind in data.index]
            metadata.index = [str(int(ind)) for ind in metadata.index]

        return data, metadata

    def graph(self, style = 'line', date = dt.date(year = 2018, month = 4, day = 1), **kwargs):
        """
        :param style: Either 'bar' or 'line'. Note that line graphs are not supported for Realtor data.
        :param date: If doing a bar graph of Zillow data, the date to graph.
        :param plot: If false, do not actually show the plot of the data.
        :param kwargs: kwargs to pass to the process_property_data function.
        :return: None
        """

        if self.geo_filter_values != 'all':
            raise TypeError('Cannot graph when data is not filtered (i.e. when self.geo_filter_values = "all"')

        data, metadata = self.process_property_data(**kwargs)

        xtick_rotation = 90
        xtick_size = 8

        if self.source == 'Zillow' and style == 'line':
            for location in self.geo_filter_values:
                ids = [region for region in data.columns if metadata.loc[region, self.geo_filter] == location]
                filtered_data = data[ids]
                filtered_data.plot(title = location + ' ,' + self.name)
                plt.show()

        elif self.source == 'Zillow':
            for location in self.geo_filter_values:
                ids = [region for region in data.columns if metadata.loc[region, self.geo_filter] == location]
                filtered_data = data.loc[date, ids].sort_values()
                filtered_data.plot(kind = 'bar', color = 'Green', title = location + ', ' + self.name, legend = False)
                plt.xticks(rotation=xtick_rotation, size = xtick_size)
                plt.show()

        elif self.source == 'Realtor':
            for location in self.geo_filter_values:
                ids = [region for region in data.index if metadata.loc[region, self.geo_filter] == location]
                filtered_data = data.loc[ids].sort_values(self.name)
                filtered_data.plot(kind = 'bar', title = location)
                plt.xticks(rotation=xtick_rotation, size = xtick_size)
                plt.show()

# Zillow inputs --
sfhomes_nbhd = Property_Input(path = 'data/Zillow Data/sfhomes_neighborhood.csv',
                        source = 'Zillow')

inventoryraw = Property_Input(path = "data/Zillow Data/inventoryraw_zip.csv",
                            source = 'Zillow')

medpricecuts = Property_Input(path = "data/Zillow Data/medpricecut_zip.csv",
                            source = 'Zillow')

allhomeprices = Property_Input(path = 'data/Zillow Data/allhomeprice_zip.csv',
                             source = 'Zillow')

# Realtor inputs --
realtor_hotness_data = Property_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                                      source = 'Realtor')

realtor_core_inventory = Property_Input(path = 'data/RDC_InventoryCoreMetrics_Zip.csv',
                                 source = 'Realtor', geo_filter_values = 'all')

realtor_core_inventory_sf = Property_Input(path = "data/RDC_InventoryCoreMetrics_Zip_sfh.csv",
                                 source = 'Realtor', geo_filter_values = 'all')

realtor_core_inventory_mf = Property_Input(path = "data/RDC_InventoryCoreMetrics_Zip_cth.csv",
                                 source = 'Realtor', geo_filter_values = 'all')
