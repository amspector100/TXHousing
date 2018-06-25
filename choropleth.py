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
dallas_regulations_path = 'data/dallas zoning standards.csv'
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
                          start_color = 'blue', end_color = 'red', show = False):

    # Get rid of nas
    gdf = gdf.loc[gdf[factor].notnull()]

    # Create colormap with caption
    min_data = gdf[factor].min()
    max_data = gdf[factor].max()

    colormap =  cm.LinearColormap(colors = [start_color, end_color], vmin = min_data, vmax = max_data)
    if scale_name is None:
        colormap.caption = layer_name
    else:
        colormap.caption = scale_name

    # Create gjson
    gjson = folium.GeoJson(
            gdf,
            show = show,
            name = layer_name,
            style_function = lambda feature: {
                'fillColor': colormap(feature['properties'][factor]),
                'color': colormap(feature['properties'][factor]),
                'weight': weight,
                'alpha': alpha,
            }
        )

    return gjson, colormap

# Calculate non-categorical values for all zip codes ----------------------------------------------------------------
def create_regulatory_layers(zoning_input, zips, regulations_path, regulation_types, broaden = True):

    # Get zoning data
    zoning_data = helpers.process_zoning_shapefile(zoning_input, broaden = broaden)

    # Process data by only considering the polygons - might fix this later, it only excludes 40/21.6K zones though.
    zoning_data = sf.process_geometry(zoning_data)

    # Get zip intersections and zip boundaries
    zoning_data = sf.zip_intersect(zoning_data, zips)

    # Calculate areas
    zoning_data['area'] = zoning_data['geometry'].area

    # Get zip geodata and fill
    zip_geodata = helpers.get_zip_boundaries()
    zip_geodata = zip_geodata.loc[zip_geodata.index.isin([str(s) for s in zips])]
    zip_geodata = zip_geodata[['geometry']]
    zip_geodata['area'] = zip_geodata['geometry'].area

    # This is probably a bit slow... oh well.
    print('Calculating features by zip code')
    for feature in tqdm(regulation_types):

        # This can definitely be optimized by calculating base zones just once instead of doing it many many times
        zoning_data[feature] = helpers.get_regulation_data(zones = zoning_data[zoning_input.feature],
                                                           regulations_path = regulations_path,
                                                           regulation_feature = feature,
                                                           fill = regulation_types[feature]) # Work on this, the fill might need to be different
        zip_feature_data = pd.Series(index = zip_geodata.index)

        for code in zip_geodata.index:

            # Get weighted average of regulation feature

            zip_feature_data[code] = zoning_data[code].multiply(zoning_data[feature]).multiply(zoning_data['area']).sum()/zip_geodata.loc[code, 'area']

        zip_geodata[feature] = zip_feature_data.copy()

    return zip_geodata

# Add demand noncategorical data to zip geodata
def add_demand_data(zip_geodata, demand_input, city, feature_name = None):
    """
    :param zip_geodata: Zip geodata, should already have all the geometries for all of the zip codes in the city.
    :param demand_input: A class demand_input.
    :param city: The City, i.e. 'Austin, TX' or 'Dallas, TX'
    :param feature_name: Updated feature name in case of conflicting features from different sources (i.e. 'Avg Listing Price' is in both the sfhomes and the cth homes dataset)
    :return: zip_geodata updated with the new value
    """

    if feature_name is None:
        feature_name = demand_input.feature

    data, metadata = process_demand_data(demand_input, graph = False)
    data = data.loc[[index for index in data.index if metadata.loc[index, demand_input.geo_filter] == city]]
    zip_geodata[feature_name] = data[demand_input.feature]
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
    simplified_zoning_data = sf.process_geometry(raw_zoning_data)
    simplified_zoning_data.index = np.arange(0, len(simplified_zoning_data), 1)
    simplified_zoning_data = sf.combine_all(simplified_zoning_data, 'broad_zone', max_comb = 25, ignore_features=True)
    simplified_zoning_data = simplified_zoning_data.loc[simplified_zoning_data['broad_zone'] != 'Other']
    base_zones = FeatureGroup('Base Zoning', show = False)
    categorical_choropleth(simplified_zoning_data, 'broad_zone').add_to(base_zones)

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
    downzone_fg = FeatureGroup(name='Downzoned Locations', show = False)
    MarkerCluster(downzonepts).add_to(downzone_fg)

    constructionpts = [retrieve_coords(point) for point in construction['centroids']]
    construction_fg = FeatureGroup(name = 'New Construction', show = False)
    MarkerCluster(constructionpts).add_to(construction_fg)

    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished creating downzone and construction color markers, took {}. Now creating historic zones markers.'.format(time.time() - time0))

    # Use texas historical sites data for national zones - - - - -
    national_hd_fg = FeatureGroup(name = 'National Historic Registry Zones', show = False)
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
    local_hd_fg = FeatureGroup(name = 'Local Historic Districts', show = False)

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
    landmark_fg = FeatureGroup('Historic Landmarks', show = False)

    landmark_data = gpd.read_file(landmark_path)
    landmark_locations = [retrieve_coords(point) for point in landmark_data['geometry']]
    MarkerCluster(landmark_locations).add_to(landmark_fg)

    # ---------------------------------------------------See print statement-------------------------------------------
    print('Finished creating historic layers, took {}. Now creating demand layers.'.format(time.time() - time0))

    zip_geodata = create_regulatory_layers(zoning_input, austin_zips, austin_regulations_path, austin_regulation_types, broaden = False)

    # Get demand data
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_med_dom_cbsa, city = 'Austin, TX')
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_avg_sf_price, city = 'Austin, TX', feature_name = 'sf_avg_listing')
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_avg_cth_price, city = 'Austin, TX', feature_name = 'mf_avg_listing')



    print('Retreiving basemap')
    basemap = folium.Map([zoning_input.long, zoning_input.lat], zoom_start=zoning_input.zoom)


    print('Adding things to basemap')
    # Add regulatory and demand factors
    for factor in zip_features_dic:
        print('here', zip_geodata[factor])
        time.sleep(1)
        gjson, colormap = continuous_choropleth(zip_geodata, factor, zip_features_dic[factor][0], zip_features_dic[factor][1], show = False)
        colormap.add_to(basemap)
        gjson.add_to(basemap)
        BindColormap(gjson, colormap).add_to(basemap)

    # Add to basemap
    base_zones.add_to(basemap)
    construction_fg.add_to(basemap)
    downzone_fg.add_to(basemap)
    national_hd_fg.add_to(basemap)
    local_hd_fg.add_to(basemap)
    landmark_fg.add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save('Figures/Bucket 2/Austin_Mastermap.html')

 # Final graph for Dallas ---------------------------------------------------------------------------------------
def dallas_final_graph(zoning_input, zip_features_dic, included_counties = ['Dallas']):
    """
    Final graph for dallas.
    :param zoning_inputs:
    :param zip_features_dic:
    :param included_counties: List of counties, i.e. ['Dallas', 'Denton']. Will only graph base zones in those counties.
    :return:
    """

    time0 = time.time()

    # Reading initial data

    # Base zones JUST for Dallas municipality. This makes it easier to begin with.
    dallas_zoning_data = process_zoning_shapefile(dallas_inputs, broaden = True)
    base_zones = FeatureGroup('Base Zoning', show = False)
    categorical_choropleth(dallas_zoning_data, factor = 'broad_zone').add_to(base_zones)

    # Base zones for all of north texas, simplified

    # Step 1: Base zoning. Only consider in some counties and then simplify to save memory and reduce HTML file size.
    #print('In Dallas final graph call, starting to work on base zoning at time {}'.format(time.time() - time0))

    #raw_zoning_data = process_zoning_shapefile(zoning_input, broaden = True)
    #raw_zoning_data = raw_zoning_data.loc[raw_zoning_data['COUNTY'].isin(included_counties)]
    #simplified_zoning_data = sf.process_geometry(raw_zoning_data)
    #simplified_zoning_data = simplified_zoning_data.loc[simplified_zoning_data['broad_zone'].isin(['Single Family', 'Multi-Family'])]
    #simplified_zoning_data = sf.process_geometry(raw_zoning_data)
    #simplified_zoning_data.index = np.arange(0, len(simplified_zoning_data), 1)
    #simplified_zoning_data = sf.combine_all(simplified_zoning_data, 'broad_zone', max_comb = 25, ignore_features=True)
    #simplified_zoning_data = simplified_zoning_data.loc[simplified_zoning_data['broad_zone'] != 'Other']
    #base_zones = FeatureGroup('Base Zoning')
    #categorical_choropleth(simplified_zoning_data, 'broad_zone').add_to(base_zones)

    # Step 2: Historic subdistricts
    print('In Dallas final graph call, starting to work on historic districts at time {}'.format(time.time() - time0))

    # Start with conservation districts

    conservation_fg = FeatureGroup("Conservation Districts", show = False)

    cd_districts = dallas_zoning_data.loc[dallas_zoning_data['LONG_ZONE_'].apply(lambda x: x[0:2]) == 'CD']
    cd_districts = cd_districts.to_crs({'init': 'epsg:4326'})

    folium.GeoJson(
        cd_districts,
        style_function=lambda feature: {
            'fillColor': 'Blue',
            'color': 'Blue',
            'weight': global_weight,
            'fillOpacity': global_alpha,
        }
    ).add_to(conservation_fg)

    # Next do historic overlays
    historic_overlay_fg = FeatureGroup('Historic Overlays', show = False)
    historic_overlay_path = "data/Zoning Shapefiles/HistoricOverlay_Dallas/HistoricOverlay.shp"
    hist_overlays = gpd.read_file(historic_overlay_path)
    hist_overlays = hist_overlays.to_crs({'init': 'epsg:4326'})

    folium.GeoJson(
        hist_overlays,
        style_function=lambda feature: {
            'fillColor': 'Purple',
            'color': 'Purple',
            'weight': global_weight,
            'fillOpacity': global_alpha,
        }
    ).add_to(historic_overlay_fg)

    # Now do historic subdistricts
    historic_subdistricts_fg = FeatureGroup('Historic Subdistricts', show = False)
    historic_subdistricts_path = "data/Zoning Shapefiles/HistoricSubdistricts_Dallas/HistoricSubdistricts.shp"
    historic_subdistricts = gpd.read_file(historic_subdistricts_path)
    historic_subdistricts = historic_subdistricts.to_crs({'init': 'epsg:4326'})

    folium.GeoJson(
        historic_subdistricts,
        style_function=lambda feature: {
            'fillColor': 'Red',
            'color': 'Red',
            'weight': global_weight,
            'fillOpacity': global_alpha,
        }
    ).add_to(historic_subdistricts_fg)

    # National historic zones
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

    # Now create regulatory layers -----------------------------------------------------------------------------------
    zip_geodata = create_regulatory_layers(dallas_inputs, dallas_zips, dallas_regulations_path, dallas_regulation_types, broaden = False)

    # Get data, subset to Dallas, get geodata

    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_med_dom_cbsa, city = 'Dallas, TX')
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_avg_sf_price, city = 'Dallas, TX', feature_name = 'sf_avg_listing')
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_avg_cth_price, city = 'Dallas, TX', feature_name = 'mf_avg_listing')


    # Final list of featuregroups: base_zones, conservation_fg, historic_overlay_fg
    print('Retreiving Dallas basemap')
    basemap = folium.Map([zoning_input.long, zoning_input.lat], zoom_start=zoning_input.zoom)


    print('Adding things to Dallas basemap')
    # Add regulatory and demand factors
    for factor in zip_features_dic:
        gjson, colormap = continuous_choropleth(zip_geodata, factor, zip_features_dic[factor][0], zip_features_dic[factor][1], show = False)
        colormap.add_to(basemap)
        gjson.add_to(basemap)
        BindColormap(gjson, colormap).add_to(basemap)

    # Add other featuregroups
    base_zones.add_to(basemap)
    conservation_fg.add_to(basemap)
    historic_overlay_fg.add_to(basemap)
    historic_subdistricts_fg.add_to(basemap)
    LayerControl().add_to(basemap)

    basemap.save('Figures/Bucket 2/Dallas_Mastermap.html')


nbhd_boundaries_path = "data\Zillow Data\ZillowNeighborhoods-TX\ZillowNeighborhoods-TX.shp"
nbhd_boundaries = gpd.read_file(nbhd_boundaries_path)
nbhd_boundaries = nbhd_boundaries.loc[nbhd_boundaries['County'] == 'Dallas'].set_index('Name')

zip_geodata = helpers.get_zip_boundaries()
zip_geodata = zip_geodata.loc[dallas_zips]



zip_geodata.plot(color = 'red', alpha = 0.5)
nbhd_boundaries.plot(color = 'blue', alpha = 0.5)
plt.show()


#dallas_final_graph(north_texas_inputs, dallas_zip_features_dic, included_counties = ['Dallas'])







