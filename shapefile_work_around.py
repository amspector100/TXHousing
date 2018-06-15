import numpy as np
import matplotlib.pyplot as plt
import shapely
import json
import pandas as pd
import geopandas as gdp
import shapefile as shp
import os
import time
import folium
from shapely.geometry.multilinestring import *
import copy

# Just format shp_inputs quickly
class zoning_inputs():

    def __init__(self, path, feature, separator, proj4string, base_zones, lat = 0, long = 0, zoom = 10,
                 title = '', xlims = [None, None], ylims = [None, None]):
        self.path = path
        self.feature = feature
        self.separator = separator
        self.proj4string = proj4string
        self.base_zones = base_zones
        self.lat = lat
        self.long = long
        self.zoom = zoom
        self.title = title
        self.xlims = xlims
        self.ylims = ylims

# Globals
# Downloaded from http://data-nctcoggis.opendata.arcgis.com/datasets/2015-land-use
north_texas_inputs = zoning_inputs(path = '2015_North_Texas_Land_Use/2015_Land_Use.shp',
                                   feature = 'CATEGORY',
                                   separator = 'no_separator',
                                   proj4string = 'EPSG:4326',
                                   base_zones = { 'Single Family':['Single family'],
                                               'Multi-Family':['Mixed use', 'Multi-family'],
                                               'Other Residential':['Residential acreage'],
                                               'Other':[]},
                                   long = 32.7767,
                                   lat=-96.7970,
                                   zoom = 10,
                                   title = 'Base Zones in and around Dallas, Texas')
                                   #xlims = [-97.3, -96.1],
                                   #ylims = [32.25, 33.25])
# Downloaded from https://data-nctcoggis.opendata.arcgis.com/datasets/81916863ab394786ab5caaa731f5ac36_4?geometry=-104.084%2C30.991%2C-83.342%2C34.229
nt_highways = zoning_inputs(path = "NT_Highways_2017\Highways_2017.shp",
                            feature = 'unsure',
                            separator = 'no_separator',
                            proj4string = 'EPSG:4326',
                            base_zones = {},
                            long = north_texas_inputs.long,
                            lat = north_texas_inputs.lat,
                            zoom = north_texas_inputs.zoom)

# Downloaded from https://data-nctcoggis.opendata.arcgis.com/datasets/counties-
nt_counties = zoning_inputs(path = "nt_counties\Counties_.shp",
                            feature = '',
                            separator = '',
                            proj4string = 'EPSG:4326',
                            base_zones = {},
                            long = north_texas_inputs.long,
                            lat = north_texas_inputs.lat,
                            zoom = north_texas_inputs.zoom)



# Downloaded from https://data.austintexas.gov/Locations-and-Maps/Zoning/5rzy-nm5e
austin_inputs = zoning_inputs(path = "austin_zoning/geo_export_571668ee-52f1-4ac9-a4e0-3b8bb348eae7.shp",
                              feature = 'zoning_zty',
                              separator = '-',
                              proj4string = 'EPSG:4326',
                              base_zones = {1:['SF'], 2:['MF'],
                                            3:['MH', 'RR', 'LA'], 4:[]},
                              long = 30.267,
                              lat = -97.743,
                              zoom = 9,
                              title = 'Base Zones in Austin, Texas')

os.chdir("C:/Users/amspe/Documents/R/MI2018/TXHousing/data/Zoning Shapefiles")


def simple_process_shapefile(input, plot = True,
                             highway_inputs = nt_highways, plot_highways = True,
                             boundary_inputs = nt_counties, plot_boundaries = True):
    """"
    Processes shapefile. Everything should be clear except conversions, which should be a dictionary
    which maps base zone types (i.e. sf, mf) to a list of strings which count as that base zone type. See
    austin_inputs['conversions'] for an example.
    """

    # Begin timing
    time0 = time.time()

    # Make sure input is correct
    if isinstance(input, zoning_inputs) == False:
        warning('Error, Input must be of class zoning_inputs')
        return None

    # Read data and process zone codes
    print('Reading file')
    raw_data = gdp.read_file(input.path)
    print('Finished reading file, took {}'.format(time.time() - time0))

    def get_zone(text):
        split_text = text.split(input.separator)[0]
        for key in input.base_zones:
            if split_text in input.base_zones[key]:
                return(key)
        return([key for key in input.base_zones][-1])

    raw_data['zone_code'] = raw_data[input.feature].apply(get_zone)
    print('Finished processing zones, took {}'.format(time.time() - time0))

    # In the future, you could easily plot basemap shapefiles underneath this

    # Plot
    if plot:

        # Make color dictionary
        color_list = plt.cm.tab10(np.linspace(0, 1, len(input.base_zones)))
        color_dic = {}
        for color, key in zip(color_list, [key for key in input.base_zones]):
            color_dic[key] = color

        # Plot objects
        fig, ax = plt.subplots()
        legend_handlers = []

        # Boundaries
        # Zoning data
        print('Beginning to graph zoning data')
        for zone in [key for key in input.base_zones]:
            filtered_data = raw_data.loc[raw_data['zone_code'] == zone]
            gdp.plotting.plot_polygon_collection(ax, filtered_data['geometry'].copy(), color = color_dic[zone], alpha = 0.7, label = zone)
            legend_handlers.append(plt.scatter([], [], color=color_dic[zone]))
        print('Finished graphing zoning data, took {}'.format(time.time() - time0))

        if plot_boundaries:
            print('Reading and plotting boundary data')
            boundary_data = gdp.read_file(boundary_inputs.path)
            gdp.plotting.plot_linestring_collection(ax, boundary_data['geometry'].boundary, color = 'black', linewidths = (0.5,))
            for point, countyname in zip(boundary_data['geometry'].centroid, boundary_data['COUNTY']):
                ax.text(point.coords[:][0][0], point.coords[:][0][1], countyname, size = 6, ha = 'center')

            print('Finished handling boundary data, time is {}'.format(time.time() - time0))


        # Highways
        if plot_highways:
            print('Reading and plotting highway data')
            highway_data = gdp.read_file(highway_inputs.path)

            primary_highways = highway_data.loc[highway_data['CLASS'] == 'Primary Highway']
            print(primary_highways)
            gdp.plotting.plot_linestring_collection(ax, primary_highways['geometry'], alpha = 1,
                                                    color = 'black', linewidths = (0.3, ))

            # Plot highwaynames
            counter = {}
            for num in primary_highways['HWY_NUM'].unique():
                counter[num] = 0

            for point, highwayname in zip(primary_highways['geometry'].centroid, primary_highways['HWY_NUM']):
                if counter[num] == 0:
                    ax.text(point.coords[:][0][0], point.coords[:][0][1], highwayname, size = 7, ha = 'left')
                    counter[num] += 1
                elif counter[num] == 150:
                    counter[num] = 0
                else:
                    counter[num] += 1
                    continue

            secondary_highways = highway_data.loc[highway_data['CLASS'] == 'Secondary Highway']
            gdp.plotting.plot_linestring_collection(ax, secondary_highways['geometry'], alpha = 1,
                                                    color = 'black', linewidths = (0.2, ))
            print('Finished handling highway data, time is {}'.format(time.time() - time0))

        # Labels and limits
        ax.set_xlabel('Longitude')
        ax.set_xlim(left = input.xlims[0], right = input.xlims[1])
        ax.set_ylabel('Latitude')
        ax.set_ylim(bottom = input.ylims[0], top = input.ylims[1])
        ax.set_title(input.title)
        ax.legend(tuple(legend_handlers), tuple([key for key in input.base_zones]))


        print('Finished graphing, took {}. Beginning to show.'.format(time.time() - time0))
        plt.show()

    return raw_data

# More complicated graphing with osm basefile using folium
def layered_graphing(input, **kwargs):

    try:
        processed_data = kwargs['processed_data']
    except:
        processed_data = simple_process_shapefile(input, plot = False)

    processed_data = processed_data.to_crs({'init': 'epsg:4326'})

    # Get basemap
    print('Retreiving basemap')
    basemap = folium.Map([input.long, input.lat], zoom_start=input.zoom)

    # Plot choropleth
    time0 = time.time()
    print('Plotting choropleth')
    basemap.choropleth(geo_data = processed_data,
                       data = processed_data,
                       columns = ['shape_area', 'zone_code'],
                       key_on = 'feature.properties.shape_area',
                       fill_color = 'PuBuGn',
                       legend_name = 'Base Zone',
                       fill_opacity = 0.6,
                       line_opacity = 0)

    print('Finished with choropleth call, took {}'.format(time.time() - time0))
    return basemap

if __name__ == '__main__':

    processed_data = simple_process_shapefile(austin_inputs, plot = False, plot_boundaries=False, plot_highways=False)

#else:
    #data = simple_process_shapefile(north_texas_inputs, plot = True, title = 'Base Zoning in and around Dallas, TX')
    #small_data = data.iloc[1:1000]