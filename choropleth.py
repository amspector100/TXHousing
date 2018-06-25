# Goal is to create layered folium graphs which allow better understanding of zoning in Austin, Dallas, and (maybe) Houston
import os
import sys
import time
import copy
import warnings
from tqdm import tqdm

import datetime as dt
import numpy as np
import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import folium
from folium import FeatureGroup, LayerControl
from folium.plugins import MarkerCluster
import branca.colormap as cm

from inputs import *
import helpers
from helpers import *
import spatial_functions as sf
import nonspatialgraphing as nsg
from BindColorMap import *

# To do
# 1. Add legend functionality for categorical choropleth, somehow.

# Globals
os.chdir("C:/Users/amspe/Documents/R/MI2018/TXHousing")
austin_regulations_path = "data/austin zoning standards.csv"
global_weight = 1
global_alpha = 0.7

# Very basic choropleth functions -------------------------------------------------------------------------------------
def categorical_choropleth(gdf, factor, colors = None, quietly = False, weight = global_weight, alpha = global_alpha):
    """
    Creates categorical choropleth using Blues spectrum
    :param gdf: A geopandas geodataframe.
    :param factor: The feature you want to plot (should be categorical).
    :param colors: Colors to use in the categorical plot. Will generate colors using the tab10 colormap.
    :param quietly: If true, will not print anything.
    :return: A folium geojson.
    """

    values = gdf[factor].unique()

    # Get colors
    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(values)))
        colors = [helpers.convert_to_hex(color) for color in colors]
    elif len(colors) < len(values):
        raise IndexError('In categorical_choropleth call, the "colors" input has fewer colors than the data has unique values')

    # Get colordic and apply to data
    colordic = {value: color for value, color in zip(values, colors)}
    gdf['color'] = gdf[factor].map(colordic)
    if not quietly:
        print('Legend functionality is not available in categorical choropleths yet, so instead we print the colordic')
        print('Here it is: {}'.format(colordic))

    # Transform data, as always
    gdf = gdf.to_crs({'init': 'epsg:4326'})

    gjson = folium.GeoJson(
            gdf,
            style_function=lambda feature: {
                'fillColor': feature['properties']['color'],
                'color': feature['properties']['color'],
                'weight': weight,
                'fillOpacity': alpha,
            }
        )

    return gjson

#  continuous colormap
def continuous_choropleth(gdf, factor, layer_name, scale_name = None, weight = global_weight, alpha = global_alpha,
                          start_color = 'blue', end_color = 'red'):

    # Get rid of nas
    gdf = gdf.loc[gdf[factor].notnull()]

    # Create colormap with caption
    min_data = gdf[factor].min()
    max_data = gdf[factor].max()
    print(min_data, max_data)

    colormap =  cm.LinearColormap(colors = [start_color, end_color], vmin = min_data, vmax = max_data)
    print(min_data, max_data, colormap((max_data - min_data)/2), factor)
    if scale_name is None:
        colormap.caption = layer_name
    else:
        colormap.caption = scale_name

    for ind, row in gdf.iterrows():
        print(ind, row[factor])
        print(colormap(row[factor]))

    # Create gjson
    gjson = folium.GeoJson(
            gdf,
            name = layer_name,
            style_function = lambda feature: {
                'fillColor': colormap(feature['properties'][factor]),
                'color': colormap(feature['properties'][factor]),
                'weight': weight,
                'alpha': alpha,
            }
        )

    return gjson, colormap

# Calculate non-categorical values for all zip codes
def create_regulatory_layers(zoning_input, zips, regulations_path, broaden = True):

    # Get zoning data
    zoning_data = helpers.process_zoning_shapefile(zoning_input, broaden = broaden)

    # Process data by only considering the polygons - might fix this later, it only excludes 40/21.6K zones though.
    zoning_data = sf.process_geometry(zoning_data)

    # Get zip intersections and zip boundaries
    zoning_data = sf.zip_intersect(zoning_data, zips)

    # Get zip geodata and fill
    zip_geodata = helpers.get_zip_boundaries()
    zip_geodata = zip_geodata.loc[zip_geodata.index.isin([str(s) for s in zips])]
    zip_geodata = zip_geodata[['geometry']]

    # This is probably a bit slow... oh well.
    print('Calculating features by zip code')
    for feature in tqdm(austin_regulation_types):

        # This can definitely be optimized by calculating base zones just once instead of doing it many many times
        zoning_data[feature] = helpers.get_regulation_data(zones = zoning_data[zoning_input.feature],
                                                           regulations_path = regulations_path,
                                                           regulation_feature = feature,
                                                           fill = austin_regulation_types[feature]) # Work on this, the fill might need to be different

        zip_feature_data = pd.Series(index = zip_geodata.index)

        for code in zip_geodata.index:

            # Get weighted average of regulation feature

            zip_feature_data[code] = zoning_data[code].multiply(zoning_data[feature]).multiply(zoning_data['geometry'].area).sum()/zip_geodata.loc[code, 'geometry'].area

        zip_geodata[feature] = zip_feature_data.copy()

    return zip_geodata

# Actually getting the final objects -------------------------------AUSTIN----------------------------------------------

def final_austin_graph(zoning_input, zip_features_dic):
    """
    :param zoning_input: Of class zoning input. Will basically only work for the austin input.
    :param zip_features_dic: zip_features_dic: a dictionary of features that will end up in the zip_geodata
    geodataframe mapped to a list with two elements, the first being the layer name and the second being the
    name for the caption.
    I.e. {'meddom':['Median Days on Market', 'Median Days on Market, Data from Realtor']}
    :return:
    """


    time0 = time.time()

    # Get initial parcel data, zip data, and zoning data
    print('Reading initial data')

    # Make sure input is correct
    if isinstance(zoning_input, zoning_inputs) == False:
        warning('Error, Input must be of class zoning_inputs')
        return None

    raw_zoning_data = process_zoning_shapefile(zoning_input, broaden = True)

    # Simplify categorically for plotting
    #simplified_zoning_data = sf.process_geometry(raw_zoning_data)
    #simplified_zoning_data.index = np.arange(0, len(simplified_zoning_data), 1)
    #simplified_zoning_data = sf.combine_all(simplified_zoning_data, 'broad_zone', max_comb = 25, ignore_features=True)
    #simplified_zoning_data = simplified_zoning_data.loc[simplified_zoning_data['broad_zone'] != 'Other']
    #base_zones = FeatureGroup('Base Zoning')
    #categorical_choropleth(simplified_zoning_data, 'broad_zone').add_to(base_zones)

    parcel_data = gpd.read_file(nsg.austin_processed_parcel_data_path)
    parcel_data['centroids'] = parcel_data['geometry'].centroid
    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished reading processed parcel data and zoning data, took {}. Starting to create downzone and construction color markers.'.format(time.time() - time0))

    downzone, construction = nsg.austin_development_histogram(parcel_data, plot = False)

    # Make marker clusters
    def retrieve_coords(point):
        result = list(point.coords[:][0][0:])
        result.reverse()
        return result

    downzonepts = [retrieve_coords(point) for point in downzone['centroids']]
    downzone_fg = FeatureGroup(name='Downzoned Locations')
    MarkerCluster(downzonepts).add_to(downzone_fg)

    constructionpts = [retrieve_coords(point) for point in construction['centroids']]
    construction_fg = FeatureGroup(name = 'New Construction')
    MarkerCluster(constructionpts).add_to(construction_fg)

    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished creating downzone and construction color markers, took {}. Now creating historic zones markers.'.format(time.time() - time0))

    # Use texas historical sites data for national zones - - - - -
    national_hd_fg = FeatureGroup(name = 'National Historic Registry Zones')
    tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
    tx_hd_data = gpd.read_file(tx_hd_path)

    folium.GeoJson(
        tx_hd_data,
        style_function=lambda feature: {
            'fillColor': 'Green',
            'color': 'Green',
            'weight': global_weight,
            'fillOpacity': global_alpha,
        }
    ).add_to(national_hd_fg)

    # Get local districts - - - - - - - - -

    signature = '-HD'
    local_hd_fg = FeatureGroup(name = 'Local Historic Districts')

    # Only consider historic districts
    print('Filtering to only include historic districts')
    processed_data = raw_zoning_data.loc[[signature in text for text in raw_zoning_data[zoning_input.feature]]]
    processed_data = processed_data.to_crs({'init': 'epsg:4326'})

    folium.GeoJson(
        processed_data,
        style_function=lambda feature: {
            'fillColor': 'Blue',
            'color': 'Blue',
            'weight': global_weight,
            'fillOpacity': global_alpha,
        }
    ).add_to(local_hd_fg)

    # Get local landmarks - - - - -
    landmark_path = "data/Zoning Shapefiles/Historical Landmarks/geo_export_ce453b58-d8ca-47d4-8934-76d636d24ca3.shp"
    landmark_fg = FeatureGroup('Historic Landmarks')

    landmark_data = gpd.read_file(landmark_path)
    landmark_locations = [retrieve_coords(point) for point in landmark_data['geometry']]
    MarkerCluster(landmark_locations).add_to(landmark_fg)

    # ---------------------------------------------------See print statement-------------------------------------------
    print('Finished creating historic layers, took {}. Now creating demand layers.'.format(time.time() - time0))

    zip_geodata = create_regulatory_layers(zoning_input, austin_zips, austin_regulations_path, broaden = False)

    # Get data, subset to austin, get geodata
    meddom_data, meddom_metadata = process_demand_data(realtor_med_dom_cbsa, graph = False)
    meddom_data = meddom_data.loc[[index for index in meddom_data.index if meddom_metadata.loc[index, realtor_med_dom_cbsa.geo_filter] == 'Austin, TX']]
    zip_geodata['meddom'] = meddom_data[realtor_med_dom_cbsa.feature]

    print('Retreiving basemap')
    basemap = folium.Map([test_input.long, test_input.lat], zoom_start=test_input.zoom)


    print('Adding things to basemap')
    # Add regulatory and demand factors
    for factor in zip_features_dic:
        gjson, colormap = continuous_choropleth(zip_geodata, factor, zip_features_dic[factor][0], zip_features_dic[factor][1])
        colormap.add_to(basemap)
        gjson.add_to(basemap)
        BindColormap(gjson, colormap).add_to(basemap)




    # Add to basemap
    #base_zones.add_to(basemap)
    construction_fg.add_to(basemap)
    downzone_fg.add_to(basemap)
    national_hd_fg.add_to(basemap)
    local_hd_fg.add_to(basemap)
    landmark_fg.add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save('Figures/Bucket 2/Austin_Mastermap.html')


# Testing

test_input = austin_inputs

zip_features_dic = {'meddom':['Median Days on Market', 'Median Days on Market, Data from Realtor'],
                    'far':['Floor to Area Ratio', 'Floor to Area Ratio, Data from City of Austin'],
                    'min_lot':['Minimum Lot Size', 'Minimum Lot Size (Square Feet)'],
                    'max_build_cov':['Maximum Building Coverage', 'Maximum Building Coverage (%)']}

final_austin_graph(test_input, zip_features_dic)


#data = helpers.process_zoning_shapefile(austin_inputs, broaden = True)
#data['max_imperv'] = helpers.get_regulation_data(data[austin_inputs.feature], austin_regulations_path, 'max_imperv', fill = 100)
#print(create_regulatory_layers(austin_inputs, austin_zips, austin_regulations_path, broaden = False))



