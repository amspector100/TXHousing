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
from .BindColorMap import BindColormap



def final_houston_graph(save_path = 'Figures/Mastermaps/Houston_Mastermap.html'):
    """
    Final Houston Mastermap

    :param save_path: The path to save the html to
    :return: None
    """

    print('Retreiving Houston basemap')
    basemap = folium.Map([data_processing.zoning.houston_inputs.lat, data_processing.zoning.houston_inputs.long],
                         zoom_start=data_processing.zoning.houston_inputs.zoom)

    # Historic districts -------------------------------------------------------------------

    # National
    national_hd_fg  = FeatureGroup('National Historic Districts', show = False)
    tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
    tx_hd_data = gpd.read_file(tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston']

    folium.GeoJson(
        tx_hd_data,
        style_function=lambda feature: {
            'fillColor': 'Green',
            'color': 'Green',
            'weight': 1,
            'fillOpacity': 0.6,
        }
    ).add_to(national_hd_fg)

    national_hd_fg.add_to(basemap)

    # Local historic districts
    local_hd_fg = FeatureGroup('Local Historic Districts', show = False)
    local_hd_data = gpd.read_file(houston_historic_districts_path).to_crs({'init':'epsg:4326'})

    folium.GeoJson(
        local_hd_data,
        style_function=lambda feature: {
            'fillColor': 'Blue',
            'color': 'Blue',
            'weight': global_weight,
            'fillOpacity': global_alpha,
        }
    ).add_to(local_hd_fg)
    local_hd_fg.add_to(basemap)

    # Local historic landmarks
    local_landmarks_fg = FeatureGroup('Local Historic Landmarks', show = False)
    local_landmarks_data = gpd.read_file(houston_historic_landmarks_path).to_crs({'init':'epsg:4326'})

    make_marker_cluster(local_landmarks_data, make_centroids = False, fast = True).add_to(local_landmarks_fg)
    local_landmarks_fg.add_to(basemap)

    # Construction permit data ----------------------------------------------------------------------
    print('Processing Houston permit data')
    houston_permit_data = process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE',
                                                                   'NEW AP', 'NEW HI-'],
                                                      searchin = ['PROJ_DESC'],
                                                      earliest = 2013, latest = None)

    # Subset to only include approved permits and nonempty geometries
    houston_permit_data = houston_permit_data.loc[houston_permit_data['Approval'] == 1.0]
    houston_permit_data = sf.process_points(houston_permit_data)

    # SF construction
    sfconstruction = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW S.F.', 'NEW SF', 'NEW TOWNHOUSE', 'NEW SINGLE']))]
    sf_cons_fg = FeatureGroup("Single Family Construction Permits (Marker)", show = True)
    make_marker_cluster(sfconstruction, make_centroids = False, fast = True).add_to(sf_cons_fg)
    sf_cons_fg.add_to(basemap)
    sfconstruction_grid = sf.make_point_grid(sfconstruction)
    gjson, colormap = continuous_choropleth(sfconstruction_grid, factor = 'value',
                                              layer_name='Single Family Residential Construction (Choropleth)',
                                              scale_name = 'Number of New Single Family Home Construction Permits in Area (2013-2018)',
                                              show = True)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # Heatmap
    heatmap(sfconstruction,
                       name='Single Family Construction, 2011-2016 (Heatmap)',
                       show=False,
                       radius=13,
                       min_opacity=0.5,
                       max_val=1).add_to(basemap)

    # MF construction
    mfconstruction = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW AP', 'NEW HI-']))]
    mf_cons_fg = FeatureGroup("Multifamily Construction Permits (Marker)", show = False)
    make_marker_cluster(mfconstruction, make_centroids=False, fast=True).add_to(mf_cons_fg)
    mf_cons_fg.add_to(basemap)
    mfconstruction_grid = sf.make_point_grid(mfconstruction)
    gjson, colormap = continuous_choropleth(mfconstruction_grid, factor='value',
                                              layer_name='Multifamily Residential Construction (Choropleth)',
                                              scale_name='Number of new Multifamily Construction Permits in Area (2013-2018)',
                                            mid_color='red', end_color='#660000', show=False)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # Heatmap
    heatmap(mfconstruction, name='Multifamily Construction, 2011-2016 (Heatmap)',
                       show=False,
                       radius=13,
                       min_opacity=0.5).add_to(basemap)

    # Add dark layer for visualization, layer control, then save
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save()