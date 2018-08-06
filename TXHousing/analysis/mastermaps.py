import os
import sys
import time
import copy
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd

import folium
from folium import FeatureGroup, LayerControl

# Import other submodules
from .. import utilities
from .. import data_processing
from . import choropleth

# Paths for historic districts - thankfully these mostly require no real processing
tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
houston_historic_districts_path = "data/Houston_Historic_Protections/HISTORIC_DISTRICTS_CITY.shp"
houston_historic_landmarks_path = "data/Houston_Historic_Protections/HISTORICAL_SITES.shp"
dallas_historic_overlay_path = "data/Zoning Shapefiles/HistoricOverlay_Dallas/HistoricOverlay.shp"
dallas_historic_subdistricts_path = "data/Zoning Shapefiles/HistoricSubdistricts_Dallas/HistoricSubdistricts.shp"


def create_austin_mastermap(save_path = 'Figures/Mastermaps/Austin_Mastermap_v2.html'):

    time0 = time.time()

    basemap = folium.Map([data_processing.zoning.austin_inputs.lat, data_processing.zoning.austin_inputs.long],
                         zoom_start=data_processing.zoning.austin_inputs.zoom)

    # Single family construction - - - - - - - - -

    # Marker
    sf_cons_fg = FeatureGroup('Single Family Residential Construction (Markers)', show = True)
    sfconstruction = process_austin_permit_data(searchfor = ['101 single family houses'], earliest = 2013, permittypedesc = 'Building Permit', workclass = 'New')
    make_marker_cluster(sfconstruction, make_centroids=False, fast=True).add_to(sf_cons_fg)
    sf_cons_fg.add_to(basemap)

    # Choropleth
    sfconstruction_grid = sf.make_point_grid(sfconstruction)
    gjson, colormap = continuous_choropleth(sfconstruction_grid, factor = 'value',
                                              layer_name='Single Family Residential Construction (Choropleth)',
                                              scale_name = 'Number of new Single Family Home Construction Permits in Area', show = True)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # Heatmap
    heatmap(sfconstruction,
                       name='Single Family Construction, 2013-2018 (Heatmap)',
                       show=False,
                       radius=13,
                       min_opacity=0.5,
                       max_val=1).add_to(basemap)


    # Multifamily construction - - - - - - - - - -

    # Marker
    mf_cons_fg = FeatureGroup('Multifamily Residential Construction (Markers)', show = False)
    mfconstruction = process_austin_permit_data(searchfor=    ['103 two family bldgs',
                                                               '104 three & four family bldgs',
                                                               '105 five or more family bldgs'],
                                                permittypedesc='Building Permit',
                                                workclass='New',
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

    # Heatmap

    heatmap(mfconstruction, name='Multifamily Construction, 2013-2018 (Heatmap)',
                       show=False,
                       radius=13,
                       min_opacity=0.5).add_to(basemap)

    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished permitting layers, took {}. Now creating historic zones markers.'.format(time.time() - time0))

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

    # Add areas with waived parking requirement - this is presently untested
    low_parking_areas = raw_zoning_data.loc[raw_zoning_data['zoning_zty'].str.contains(('|').join(['CBD', 'DMU', 'PD']))]
    low_parking_gjson = categorical_choropleth(low_parking_areas, factor = 'zoning_zty')
    low_parking_fg = FeatureGroup('Areas with Waived Parking Requirements', show = False)
    low_parking_gjson.add_to(low_parking_fg)


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
    national_hd_fg.add_to(basemap)
    local_hd_fg.add_to(basemap)
    landmark_fg.add_to(basemap)
    low_parking_fg.add_to(basemap)

    # Add dark layer for visualization and layer control
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)

    print("Saving basemap, time is {}".format(time.time() - time0))
    basemap.save('Figures/Mastermaps/Austin_Mastermap.html')




def create_dallas_mastermap(save_path = 'Figures/Mastermaps/Dallas_Mastermap_v2.html'):
    """
    Final graph for dallas.
    :param save_path: The path at which to save the final HTML.
    :return: None
    """

    time0 = time.time()

    basemap = folium.Map([data_processing.zoning.dallas_inputs.lat, data_processing.zoning.dallas_inputs.long],
                         zoom_start=data_processing.zoning.dallas_inputs.zoom)

    # Reading initial data
    print('Finished retrieving basemap, took {}. Reading initial data'.format(time.time() - time0))

    base_zones = FeatureGroup('Base Zoning', show = False)

    # Base zones just for Dallas municipality.
    regulation_features = ['height', 'inv_density', 'far']
    regulation_feature_names = ['Maximum Height', 'Minimum Lot Size (Sq Ft)', 'Maximum Floor to Area Ratio']

    dallas_zoning_data = data_processing.zoning.dallas_inputs.process_zoning_shapefile(regulation_features = regulation_features)
    choropleth.categorical_choropleth(dallas_zoning_data, factor = 'broad_zone', colors = ['#D62728', '#17BeCf', '#1f77B4', '#E377C2'], name = 'Base Zoning', basemap = basemap)
    for feature, name in zip(regulation_features, regulation_feature_names):
        print('Working with {}'.format(feature))
        choropleth.continuous_choropleth(dallas_zoning_data.loc[dallas_zoning_data[feature].notnull()],
                                         factor = feature, layer_name=name,
                                         show = False, basemap = basemap)

    choropleth.polygon_layer(dallas_zoning_data.loc[dallas_zoning_data['base_zone'] == 'CA-1'],
                             color = 'Blue',
                             name = 'Areas with Weaker Parking Requirements (CA-1 District)',
                             basemap = basemap)

    # Step 2: Historic subdistricts --
    print('In Dallas final graph call, starting to work on historic districts at time {}'.format(time.time() - time0))

    # Start with conservation districts
    cd_districts = dallas_zoning_data.loc[dallas_zoning_data['LONG_ZONE_'].apply(lambda x: x[0:2]) == 'CD']
    choropleth.polygon_layer(cd_districts, color = 'Blue', name = 'Conservation Districts', basemap = basemap)

    # Next do historic overlays
    hist_overlays = gpd.read_file(dallas_historic_overlay_path).to_crs({'init': 'epsg:4326'})
    choropleth.polygon_layer(hist_overlays, color = 'Purple', name = 'Historic Overlays', basemap = basemap)

    # Now do historic subdistricts
    historic_subdistricts = gpd.read_file(dallas_historic_subdistricts_path).to_crs({'init': 'epsg:4326'})
    choropleth.polygon_layer(historic_subdistricts, color = 'Red', name = 'Historic Subdistricts', basemap = basemap)

    # Now do national historic zones
    tx_hd_data = gpd.read_file(tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Dallas']
    choropleth.polygon_layer(tx_hd_data, name = 'National Historic Districts', color = 'Green', basemap = basemap)

    # Step 3: Property/regulatory layers -----------------------------------------------------------------------------
    print('In Dallas final graph call, starting to work on property layers at time {}'.format(time.time() - time0))

    zipdata = data_processing.boundaries.ZipBoundaries(ziplist = data_processing.boundaries.dallas_zips)
    zipdata.add_property_data(data_processing.property.realtor_hotness_data, rsuffix = '_hotness')
    zipdata.add_property_data(data_processing.property.realtor_core_inventory_sf, rsuffix = '_sf')
    zipdata.add_property_data(data_processing.property.realtor_core_inventory_mf, rsuffix = '_mf')
    for feature, name in zip(['Views Per Property (vs CBSA)', 'Hotness Rank ', 'Median Listing Price_sf', 'Median Listing Price_mf'],
                             ['Views Per Property (Mean-Centered)', 'Realtor Hotness Rank', 'Median Listing Price for Single Family Homes',
                                                                                 'Median Listing Price for Multifamily Homes']):
        choropleth.continuous_choropleth(zipdata.data, factor = feature, layer_name = name, show = False, basemap = basemap)


    # Construction
    all_construction = data_processing.permit.get_corrected_dallas_permit_data()

    # Construction -- single family - - - - - -

    # Marker
    sfconstruction = all_construction.loc[all_construction['Permit Type'] == 'Building (BU) Single Family  New Construction']
    choropleth.make_marker_cluster(sfconstruction, make_centroids = False, fast = True,
                                   name = "Single Family Construction Permits (Marker)", basemap = basemap)
    sfconstruction_grid = utilities.simple.make_point_grid(sfconstruction)
    choropleth.continuous_choropleth(sfconstruction_grid, factor = 'value',
                                      layer_name='Single Family Residential Construction (Choropleth)',
                                      scale_name = 'Number of New Single Family Home Construction Permits in Area (2011-2016)',
                                      show = False, basemap = basemap)
    choropleth.heatmap(sfconstruction, name='Single Family Construction, 2011-2016 (Heatmap)', show=False, radius=13, min_opacity=0.5, max_val=1, basemap = basemap)


    # Construction -- Multifamily - - - - - -
    mfconstruction = all_construction.loc[all_construction['Permit Type'] == 'Building (BU) Multi Family  New Construction']
    choropleth.make_marker_cluster(mfconstruction, make_centroids=False, fast=True,
                                   name = "Multifamily Construction Permits (Marker)", basemap = basemap)
    mfconstruction_grid = utilities.simple.make_point_grid(mfconstruction)
    choropleth.continuous_choropleth(mfconstruction_grid, factor='value',
                                    layer_name='Multifamily Residential Construction (Choropleth)',
                                    scale_name='Number of new Multifamily Construction Permits in Area (2011-2016)',
                                    show=False, basemap = basemap)
    # Heatmap
    choropleth.heatmap(mfconstruction, name='Multifamily Construction, 2011-2016 (Heatmap)', radius=13, min_opacity=0.5, basemap = basemap)


    # Add dark basemap
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save(save_path)


def create_houston_mastermap(save_path = 'Figures/Mastermaps/Houston_Mastermap_v2.html'):
    """
    Final Houston Mastermap

    :param save_path: The path to save the html to
    :return: None
    """

    basemap = folium.Map([data_processing.zoning.houston_inputs.lat, data_processing.zoning.houston_inputs.long],
                         zoom_start=data_processing.zoning.houston_inputs.zoom)

    # Historic districts -------------------------------------------------------------------

    # National
    tx_hd_data = gpd.read_file(tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston']
    choropleth.polygon_layer(tx_hd_data, name = 'National Historic Districts', color = 'Green', basemap = basemap)

    # Local historic districts
    local_hd_data = gpd.read_file(houston_historic_districts_path).to_crs({'init':'epsg:4326'})
    choropleth.polygon_layer(tx_hd_data, name = 'Local Historic Districts', color = 'Blue', basemap = basemap)

    # Local historic landmarks
    local_landmarks_data = gpd.read_file(houston_historic_landmarks_path).to_crs({'init':'epsg:4326'})
    choropleth.make_marker_cluster(local_landmarks_data, make_centroids = False, fast = True,
                                   name = 'Local Historic Landmarks', basemap = basemap)

    # Construction permit data ----------------------------------------------------------------------
    print('Processing Houston permit data')
    houston_permit_data = data_processing.permit.process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE',
                                                                                           'NEW AP', 'NEW HI-'],
                                                                              searchin = ['PROJ_DESC'],
                                                                              earliest = 2013, latest = None)

    # Subset to only include approved permits and nonempty geometries
    houston_permit_data = houston_permit_data.loc[houston_permit_data['Approval'] == 1.0]
    houston_permit_data = utilities.simple.process_points(houston_permit_data)

    # SF construction - marker cluster, choropleth, and heatmap
    sfconstruction = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW S.F.', 'NEW SF', 'NEW TOWNHOUSE', 'NEW SINGLE']))]
    choropleth.make_marker_cluster(sfconstruction, make_centroids = False, fast = True,
                                   name = "Single Family Construction Permits (Marker)", basemap = basemap)
    sfconstruction_grid = utilities.simple.make_point_grid(sfconstruction)
    choropleth.continuous_choropleth(sfconstruction_grid, factor = 'value',
                                      layer_name='Single Family Residential Construction (Choropleth)',
                                      scale_name = 'Number of New Single Family Home Construction Permits in Area (2013-2018)',
                                      show = True, basemap = basemap)

    choropleth.heatmap(sfconstruction, name='Single Family Construction, 2013-2018 (Heatmap)', show=False, radius=13,
                       min_opacity=0.5, max_val=1, basemap = basemap)

    # MF construction
    mfconstruction = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW AP', 'NEW HI-']))]
    choropleth.make_marker_cluster(mfconstruction, make_centroids=False, fast=True,
                                   name = "Multifamily Construction Permits (Marker)", basemap = basemap)
    mfconstruction_grid = utilities.simple.make_point_grid(mfconstruction)
    choropleth.continuous_choropleth(mfconstruction_grid, factor='value',
                                    layer_name='Multifamily Residential Construction (Choropleth)',
                                    scale_name='Number of new Multifamily Construction Permits in Area (2013-2018)',
                                    show=False, basemap = basemap)
    # Heatmap
    choropleth.heatmap(mfconstruction, name='Multifamily Construction, 2013-2018 (Heatmap)', radius=13,
                       min_opacity=0.5, basemap = basemap)

    # Add dark layer for visualization, layer control, then save
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save(save_path)

create_dallas_mastermap()