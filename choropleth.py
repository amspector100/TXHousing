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
from folium.plugins import MarkerCluster, FastMarkerCluster
import branca.colormap as cm

from inputs import *
import helpers
from helpers import *
from zipcode_scrape import *
import spatial_functions as sf
import nonspatialgraphing as nsg
from BindColorMap import *

# To do
# 1. Add legend functionality for categorical choropleth, somehow.

# Globals
os.chdir("C:/Users/amspe/Documents/R/MI2018/TXHousing")
global_weight = 1
global_alpha = 0.6

# Actually retrieve coords from shapely point
def retrieve_coords(point):
    result = list(point.coords[:][0][0:])
    result.reverse()
    return result

# Marker Cluster processing function
def make_marker_cluster(gdf, make_centroids = True, points_column = 'geometry', fast = False, **kwargs):
    """
    :param gdf: Geodataframe
    :param layer_name: Name of eventual Marker Cluster layer
    :param make_centroids:
    :param points_column:
    :return: MarkerCluster object
    """

    if make_centroids:
        gdf['centroids'] = gdf['geometry'].centroid
        points_column = 'centroids'

    points = [retrieve_coords(point) for point in gdf[points_column]]

    if fast:
        return FastMarkerCluster(points, **kwargs)
    else:
        return MarkerCluster(points, **kwargs)


# Very basic choropleth functions -------------------------------------------------------------------------------------
def categorical_choropleth(gdf, factor, colors = None, quietly = False, weight = global_weight, alpha = global_alpha,
                           geometry_column = 'geometry'):
    """
    Creates categorical choropleth using Blues spectrum
    :param gdf: A geopandas geodataframe.
    :param factor: The feature you want to plot (should be categorical).
    :param colors: Colors to use in the categorical plot. Will generate colors using the tab10 colormap.
    :param quietly: If true, will not print anything.
    :param weight: The weight in the style function.
    :param alpha: The alpha in the style function.
    :param geometry_column: The geometry column of the gdf. Defaults to 'geometry'.
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

    gdf = gdf[[factor, geometry_column, 'color']]

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
                          mid_color = '#00ccff', end_color = '#000066', show = False, geometry_column = 'geometry'):
    """
    :param gdf: Geodataframe
    :param factor: factor for analysis
    :param layer_name: Name of feature group layer
    :param scale_name: Name of scale
    :param weight: Weight
    :param alpha: Alpha of polygons
    :param start_color: I.e. white, for min data
    :param end_color: I.e. black, for max data
    :param show: Show by default on start
    :param geometry_column: 'geometry'
    :return:
    """

    # Get rid of nas
    gdf = gdf.loc[gdf[factor].notnull()]

    # Create colormap with caption
    min_data = gdf[factor].min()
    max_data = gdf[factor].max()

    colormap =  cm.LinearColormap(colors = ['white', mid_color, end_color], vmin = min_data, vmax = max_data).to_step(12, method = 'log')
    if scale_name is None:
        colormap.caption = layer_name
    else:
        colormap.caption = scale_name

    # Create gjson
    gdf = gdf[[factor, geometry_column]]
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

    print('Retreiving Austin basemap')
    basemap = folium.Map([zoning_input.lat, zoning_input.long], zoom_start=zoning_input.zoom)

    # Get initial parcel data, zip data, and zoning data
    print('Finished retrieving basemap, took {}. Reading initial data'.format(time.time() - time0))

    # Make sure input is correct
    if isinstance(zoning_input, zoning_inputs) == False:
        warning('Error, Input must be of class zoning_inputs')
        return None

    raw_zoning_data = process_zoning_shapefile(zoning_input, broaden = True)


    # Simplify categorically for plotting
    simplified_zoning_data = sf.process_geometry(raw_zoning_data)
    #simplified_zoning_data.index = np.arange(0, len(simplified_zoning_data), 1)
    #simplified_zoning_data = sf.combine_all(simplified_zoning_data, 'broad_zone', max_comb = 15, use_v2 = True, ignore_features = True)
    #simplified_zoning_data = sf.process_combined_result(simplified_zoning_data, factor = 'broad_zone')
    #simplified_zoning_data = simplified_zoning_data.loc[simplified_zoning_data['broad_zone'].isin(['Single Family', 'Multifamily'])]

    # Manually set CRS if necessary
    if simplified_zoning_data.crs is None:
        simplified_zoning_data.crs = {'init': 'epsg:4326'}

    base_zones = FeatureGroup('Base Zoning', show = False)
    #categorical_choropleth(simplified_zoning_data, 'broad_zone').add_to(base_zones)

    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished reading and processing zoning data, took {}. Starting to create downzone and construction color markers.'.format(time.time() - time0))

    downzone_fg = FeatureGroup(name='Downzoned Locations', show = False)
    downzone = gpd.read_file(austin_downzone_path)
    make_marker_cluster(downzone, make_centroids = True).add_to(downzone_fg)


    # Single family construction - - - - - - - - -

    # Marker
    sf_cons_fg = FeatureGroup('Single Family Residential Construction (Markers)', show = True)
    sfconstruction = process_austin_permit_data(searchfor = ['101 single family houses'], earliest = 2013)
    make_marker_cluster(sfconstruction, make_centroids=False, fast=True).add_to(sf_cons_fg)
    sf_cons_fg.add_to(basemap)


    sfconstruction_grid = sf.make_point_grid(sfconstruction)
    gjson, colormap = continuous_choropleth(sfconstruction_grid, factor = 'value',
                                              layer_name='Single Family Residential Construction (Choropleth)',
                                              scale_name = 'Number of new Single Family Home Construction Permits in Area', show = True)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # Multifamily construction - - - - - - - - - -

    # Marker
    mf_cons_fg = FeatureGroup('Multifamily Residential Construction (Markers)', show = False)
    mfconstruction = process_austin_permit_data(searchfor=    ['103 two family bldgs',
                                                               '104 three & four family bldgs',
                                                               '105 five or more family bldgs'],
                                                earliest=2013)
    make_marker_cluster(mfconstruction, make_centroids=False, fast=True).add_to(mf_cons_fg)
    mf_cons_fg.add_to(basemap)

    # Choropleth
    mfconstruction_grid = sf.make_point_grid(mfconstruction)
    gjson, colormap = continuous_choropleth(mfconstruction_grid, factor='value',
                                              layer_name='Multifamily Residential Construction (Choropleth)',
                                              scale_name='Number of new Multifamily Construction Permits in Area',
                                            mid_color='red', end_color='#660000', show=False)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished creating downzone and construction color markers, took {}. Now creating historic zones markers.'.format(time.time() - time0))

    # Use texas historical sites data for national zones - - - - -
    national_hd_fg = FeatureGroup(name = 'National Historic Registry Zones', show = False)
    tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
    tx_hd_data = gpd.read_file(tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Austin']

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


    print('Adding things to basemap')
    # Add regulatory and demand factors
    for factor in zip_features_dic:
        gjson, colormap = continuous_choropleth(zip_geodata, factor, zip_features_dic[factor][0], zip_features_dic[factor][1],
                                                mid_color='blue', end_color='red', show=False)
        colormap.add_to(basemap)
        gjson.add_to(basemap)
        BindColormap(gjson, colormap).add_to(basemap)

    # Add to basemap
    base_zones.add_to(basemap)
    downzone_fg.add_to(basemap)
    national_hd_fg.add_to(basemap)
    local_hd_fg.add_to(basemap)
    landmark_fg.add_to(basemap)

    # Add dark layer for visualization and layer control
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)

    print("Saving basemap, time is {}".format(time.time() - time0))
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

    print('Retreiving Dallas basemap')
    basemap = folium.Map([zoning_input.lat, zoning_input.long], zoom_start=zoning_input.zoom)

    # Reading initial data
    print('Finished retrieving basemap, took {}. Reading initial data'.format(time.time() - time0))

    base_zones = FeatureGroup('Base Zoning', show = False)


    # Base zones JUST for Dallas municipality. This makes it easier to begin with.
    dallas_zoning_data = process_zoning_shapefile(dallas_inputs, broaden = True)
    simplified_zoning_data = dallas_zoning_data#.loc[dallas_zoning_data['broad_zone'].isin(['Single Family', 'Multifamily'])]
    categorical_choropleth(simplified_zoning_data, factor = 'broad_zone', colors = ['#D62728', '#17BeCf', '#1f77B4', '#E377C2']).add_to(base_zones)

    # Base zones for all of north texas, simplified

    # Step 1: Base zoning. Only consider in some counties and then simplify to save memory and reduce HTML file size.
    #print('In Dallas final graph call, starting to work on base zoning at time {}'.format(time.time() - time0))

    raw_zoning_data = process_zoning_shapefile(zoning_input, broaden = True)
    raw_zoning_data = raw_zoning_data.loc[raw_zoning_data['COUNTY'].isin(included_counties)]
    simplified_zoning_data = raw_zoning_data.loc[raw_zoning_data['broad_zone'].isin(['Single Family', 'Multifamily'])]
    # --- possibly stop here

    # Or do this
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
    print('In Dallas final graph call, starting to work on regulatory and demand layers at time {}'.format(time.time() - time0))

    zip_geodata = create_regulatory_layers(dallas_inputs, dallas_zips, dallas_regulations_path, dallas_regulation_types, broaden = False)

    # Get data, subset to Dallas, get geodata

    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_med_dom_cbsa, city = 'Dallas, TX')
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_avg_sf_price, city = 'Dallas, TX', feature_name = 'sf_avg_listing')
    zip_geodata = add_demand_data(zip_geodata=zip_geodata, demand_input = realtor_avg_cth_price, city = 'Dallas, TX', feature_name = 'mf_avg_listing')


    # Downzoning and construction layers
    print('In Dallas final graph call, starting to work on downdevelopment and construction layers at time {}'.format(time.time() - time0))
    downdevelopment = gpd.read_file(dallas_downdev_path)
    downdevelopment_fg = FeatureGroup(name='Downdeveloped Locations', show = False)
    make_marker_cluster(downdevelopment, make_centroids=True).add_to(downdevelopment_fg)


    # Construction
    all_construction = get_corrected_dallas_permit_data(path = dpm_save_path)

    # Construction -- single family - - - - - -

    # Marker
    sfconstruction = all_construction.loc[all_construction['Permit Type'] == 'Building (BU) Single Family  New Construction']
    sf_cons_fg = FeatureGroup("Single Family Construction Permits (Marker)", show = True)
    make_marker_cluster(sfconstruction, make_centroids = False, fast = True).add_to(sf_cons_fg)
    sf_cons_fg.add_to(basemap)

    sfconstruction_grid = sf.make_point_grid(sfconstruction)
    gjson, colormap = continuous_choropleth(sfconstruction_grid, factor = 'value',
                                              layer_name='Single Family Residential Construction (Choropleth)',
                                              scale_name = 'Number of New Single Family Home Construction Permits in Area (2011-2016)',
                                              show = True)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # Construction -- Multifamily - - - - - -
    mfconstruction = all_construction.loc[all_construction['Permit Type'] == 'Building (BU) Multi Family  New Construction']
    mf_cons_fg = FeatureGroup('Multifamily Construction Permits (Marker)', show = False)
    make_marker_cluster(mfconstruction, make_centroids = False, fast = True).add_to(mf_cons_fg)
    mf_cons_fg.add_to(basemap)

    mfconstruction_grid = sf.make_point_grid(mfconstruction)
    gjson, colormap = continuous_choropleth(mfconstruction_grid, factor = 'value',
                                            layer_name='Multifamily Residential Construction (Choropleth)',
                                            scale_name='Number of New Multifamily Home Construction Permits in Area (2011-2016)',
                                            mid_color='red', end_color='#660000', show=False)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)


    # Final list of featuregroups: base_zones, conservation_fg, historic_overlay_fg, historic_subdistricts_fg,
    # construction_fg, downdevelopment_fg

    print('Adding things to Dallas basemap')
    # Add regulatory and demand factors
    for factor in zip_features_dic:
        gjson, colormap = continuous_choropleth(zip_geodata, factor, zip_features_dic[factor][0], zip_features_dic[factor][1],
                                                mid_color='blue', end_color='red', show=False)
        colormap.add_to(basemap)
        gjson.add_to(basemap)
        BindColormap(gjson, colormap).add_to(basemap)

    # Add other featuregroups
    base_zones.add_to(basemap)
    conservation_fg.add_to(basemap)
    historic_overlay_fg.add_to(basemap)
    historic_subdistricts_fg.add_to(basemap)
    downdevelopment_fg.add_to(basemap)

    # Add dark basemap
    folium.TileLayer('cartodbdark_matter').add_to(basemap)


    LayerControl().add_to(basemap)

    print("Saving basemap at {}".format(time.time() - time0))
    basemap.save('Figures/Bucket 2/Dallas_Mastermap.html')


#nbhd_boundaries_path = "data\Zillow Data\ZillowNeighborhoods-TX\ZillowNeighborhoods-TX.shp"
#nbhd_boundaries = gpd.read_file(nbhd_boundaries_path)
#nbhd_boundaries = nbhd_boundaries.loc[nbhd_boundaries['County'] == 'Travis'].set_index('Name')

#zip_geodata = helpers.get_zip_boundaries()
#zip_geodata = zip_geodata.loc[austin_zips]



#base = zip_geodata.plot(color = None, alpha = 0.5, edgecolor = 'red')
#nbhd_boundaries.plot(ax = base, color = None, alpha = 0.5, edgecolor = 'black')
#plt.show()

if __name__ == '__main__':

    dallas_final_graph(north_texas_inputs, dallas_zip_features_dic, included_counties = ['Dallas'])

    #final_austin_graph(austin_inputs, austin_zip_features_dic)









