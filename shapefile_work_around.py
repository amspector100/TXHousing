import numpy as np
import matplotlib.pyplot as plt
import shapely
import json
import pandas as pd
import geopandas as gpd
import shapefile as shp
import os
import time
import sys
import folium
from shapely.geometry.multilinestring import *
import copy

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
                              base_zones = {'Single Family':['SF'], 'Multifamily':['MF'],
                                            'Other Residential':['MH', 'RR', 'LA'], 'Other':[]},
                              long = 30.267,
                              lat = -97.743,
                              zoom = 9,
                              title = 'Base Zones in Austin, Texas')

austin_parcel_path = "Austin Land Database 2016/geo_export_813e97e4-7fde-4e3a-81b3-7ca9e8a89bd0.shp"
austin_parcel_feature = 'basezone'

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
    raw_data = gpd.read_file(input.path)
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
            gpd.plotting.plot_polygon_collection(ax, filtered_data['geometry'].copy(), color = color_dic[zone], alpha = 0.7, label = zone)
            legend_handlers.append(plt.scatter([], [], color=color_dic[zone]))
        print('Finished graphing zoning data, took {}'.format(time.time() - time0))

        if plot_boundaries:
            print('Reading and plotting boundary data')
            boundary_data = gpd.read_file(boundary_inputs.path)
            gpd.plotting.plot_linestring_collection(ax, boundary_data['geometry'].boundary, color = 'black', linewidths = (0.5,))
            for point, countyname in zip(boundary_data['geometry'].centroid, boundary_data['COUNTY']):
                ax.text(point.coords[:][0][0], point.coords[:][0][1], countyname, size = 6, ha = 'center')

            print('Finished handling boundary data, time is {}'.format(time.time() - time0))


        # Highways
        if plot_highways:
            print('Reading and plotting highway data')
            highway_data = gpd.read_file(highway_inputs.path)

            primary_highways = highway_data.loc[highway_data['CLASS'] == 'Primary Highway']
            print(primary_highways)
            gpd.plotting.plot_linestring_collection(ax, primary_highways['geometry'], alpha = 1,
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
            gpd.plotting.plot_linestring_collection(ax, secondary_highways['geometry'], alpha = 1,
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

    # Get basemap
    print('Retreiving basemap')
    basemap = folium.Map([input.long, input.lat], zoom_start=input.zoom)

    # Get the data
    try:
        processed_data = kwargs['processed_data']
    except:
        processed_data = simple_process_shapefile(input, plot = False)

    # Perhaps filter it
    try:
        filter = kwargs['filter']
        print('Filtering data to only include the {} base zone'.format(filter))
        processed_data = processed_data.loc[processed_data['zone_code'] == filter, ]
        plot_polygons = True

    except KeyError:
        print('Not filtering the data - creating colors for each zone')
        plot_polygons = False

        # Get colors
        try:
            colors = kwargs['colors']
        except:
            colors = plt.cm.tab10(np.linspace(0, 1, len(input.base_zones)))
            colors = [convert_to_hex(color) for color in colors]

        colordic = {base_zone:color for base_zone, color in zip([key for key in input.base_zones], colors)}
        print('Legend functionality is not available, so instead we print the colordic')
        print('Here it is: {}'.format(colordic))

        processed_data['color'] = processed_data['zone_code'].map(colordic)

    processed_data = processed_data.to_crs({'init': 'epsg:4326'})



    # Plot choropleth
    time0 = time.time()
    if plot_polygons:
        print('Plotting polygons')
        folium.GeoJson(
            processed_data,
            style_function=lambda feature: {
                'fillColor': 'Blue',
                'color': 'Blue',
                'weight': 1,
                'fillOpacity': 0.5,
            }
        ).add_to(basemap)

    else:
        print('Plotting choropleth')
        folium.GeoJson(
            processed_data,
            style_function=lambda feature: {
                'fillColor': feature['properties']['color'],
                'color': feature['properties']['color'],
                'weight': 1,
                'fillOpacity': 0.5,
            }
        ).add_to(basemap)


    print('Finished with plotting, took {}'.format(time.time() - time0))
    return basemap

# Naive attempt
def plot_historic_districts(input, signature = '-H', savename = False, plot_national_districts = True, use_nris_data = True,
                            plot_landmarks = True, landmarks_path = False, marker_cluster = True):
    """
    :param input: A "Zoning_Input" object.
    :param signature: A string used to detect whether a zoning area is a historic district.
    :return: None, but it will save a graph
    """
    
    from folium import FeatureGroup, LayerControl

    # Begin timing
    time0 = time.time()

    print('Retreiving basemap')
    basemap = folium.Map([input.long, input.lat], zoom_start=input.zoom)

    # National districts ----------------------------------------------------------

    if plot_national_districts:

        nat_dists_fg = FeatureGroup(name = 'National Registry Historic Districts')

        if use_nris_data:

            print("""You are using NRIS data, which is a bad idea, because this data has been simplified too much,
            and it makes the polygons wacky. Don't say I didn't warn you""")

            # These are not inputs, they're constant
            nat_hd_path = 'NRIS_CR_Standards_Public.gdb'
            nris_cdist = 'NRIS_crdist_py.shp'

            data = gpd.read_file(nat_hd_path, layer='NRIS_MAIN')
            data = data.loc[data['CITY'].isin(['Austin', 'Houston', 'Dallas'])]
            data = data.set_index('REFNUM')
            data = data[[column for column in data.columns if column != 'geometry']]  # Get rid of duplicate column

            # Work around fiona bug
            geodata = gpd.read_file(nris_cdist)
            geodata = geodata.set_index('NRIS_Rf')
            geodata = geodata[['geometry']]  # We don't need the other data in this layer
            national_districts = geodata.join(data, how='inner').to_crs({'init': 'epsg:4326'})

            print(national_districts)

            print('Plotting national districts')
            folium.GeoJson(
                national_districts,
                style_function=lambda feature: {
                    'fillColor': 'Red',
                    'color': 'Red',
                    'weight': 1,
                    'fillOpacity': 0.5,
                }
            ).add_to(nat_dists_fg)

        else:

            # Use texas historical sites data
            tx_hd_path = "NationalRegisterPY_shp/NationalRegisterPY.shp"
            tx_hd_data = gpd.read_file(tx_hd_path)

            print('Plotting national districts')
            folium.GeoJson(
                tx_hd_data,
                style_function=lambda feature: {
                    'fillColor': 'Green',
                    'color': 'Green',
                    'weight': 1,
                    'fillOpacity': 0.5,
                }
            ).add_to(nat_dists_fg)

        nat_dists_fg.add_to(basemap)

    # Local districts -----------------------------------------------------------------

    local_dists_fg = FeatureGroup(name = 'Local Historic Districts')

    # Make sure input is correct
    if isinstance(input, zoning_inputs) == False:
        warning('Error, Input must be of class zoning_inputs')
        return None

    # Read data and process zone codes
    print('Reading local HD file')
    raw_data = gpd.read_file(input.path)
    print('Finished reading file, took {}'.format(time.time() - time0))

    # Only consider historic districts
    print('Filtering to only include historic districts')
    raw_data[input.feature].to_csv('HDistricts.txt')
    raw_data = raw_data.loc[[signature in text for text in raw_data[input.feature]]]
    processed_data = raw_data.to_crs({'init': 'epsg:4326'})

    print('Plotting local districts')
    folium.GeoJson(
        processed_data,
        style_function=lambda feature: {
            'fillColor': 'Blue',
            'color': 'Blue',
            'weight': 1,
            'fillOpacity': 0.5,
        }
    ).add_to(local_dists_fg)

    local_dists_fg.add_to(basemap)

    # Specific historical landmarks
    if plot_landmarks:

        landmark_fg = FeatureGroup('Historic Landmarks')

        # This needs to be an input at some point
        landmark_data = gpd.read_file(landmarks_path)

        if marker_cluster:

            def retrieve_coords(point):
                result = list(point.coords[:][0][0:])
                result.reverse()
                return result

            locations = [retrieve_coords(point) for point in landmark_data['geometry']]
            print(locations)

            print('Plotting the marker cluster')
            from folium.plugins import FastMarkerCluster
            marker_cluster = FastMarkerCluster(locations).add_to(landmark_fg)

        else:

            print('Plotting landmarks')

            folium.GeoJson(landmark_data).add_to(landmark_fg)

        landmark_fg.add_to(basemap)

    LayerControl().add_to(basemap)

    if savename is False:
        basemap.save('HistoricDistricts.html')
    else:
        print('Saving under {}'.format(savename))
        basemap.save(savename)

def plot_historic_districts_parcel(input, cache = False, cached = True):

    # Now plot histogram of historic district data
    if cached:
        pass



if __name__ == '__main__':
    landmarks_path = "Historical Landmarks/geo_export_ce453b58-d8ca-47d4-8934-76d636d24ca3.shp"

    plot_historic_districts(austin_inputs, landmarks_path = landmarks_path, use_nris_data=False, signature = '-HD', savename = 'Austin_All_Hist_Districts.html')


    sys.exit()


    time0 = time.time()
    print('Starting to read Austin parcel data')
    austin_parcel = gpd.read_file(austin_parcel_path)
    print('Finished, took {} seconds'.format(time.time() - time0))
    print(austin_parcel.columns)

    for_histogram = austin_parcel[['basezone', 'market_val']]
    for_histogram.to_csv('HDHistogramData.csv')





    # Old stuff -------------------------------------------------------------------------------------------
    #plot_historic_districts(austin_inputs, 'Austin')


    #simple_process_shapefile(austin_inputs, plot = False)

    #filter = 'Multifamily'
    #mymap = layered_graphing(austin_inputs, filter = filter)
    #print('Saving the map')
    #mymap.save(filter + 'Austinmap.html')