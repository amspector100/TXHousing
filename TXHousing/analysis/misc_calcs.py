import time
import numpy as np
import pandas as pd
import geopandas as gpd
from .. import utilities
from ..data_processing import boundaries, zoning
from functools import reduce

# Special setback and minimum lot zone calculations ----------------------------- (Houston) ----------------------------------------------------- Special setback and minimum lot zone calculations
def calculate_special_zoning_areas_and_populations():
    """Calculates the area in square miles and population of (1) zones which have special minimum lot regulations and
    (2) zones which have special setback regulations. """

    special_setbacks = gpd.read_file(zoning.houston_spec_setbacks_path)
    special_lots = gpd.read_file(zoning.houston_spec_min_lot_path)

    # Find areas
    special_setbacks = utilities.measurements.get_area_in_units(special_setbacks) # Remember, this defaults to miles
    special_lots = utilities.measurements.get_area_in_units(special_lots)
    print('Special setbacks area in square miles: ', special_setbacks['area'].sum())
    print('Special lots area in square miles:', special_lots['area'].sum())

    # Find population in the particular areas - B01001e1 is the column which counts population in he block data
    block_data = boundaries.BlockBoundaries(['X01_AGE_AND_SEX'], cities = 'Houston', get_percent_residential = False)
    special_lots = block_data.push_features(special_lots, features = 'B01001e1')
    special_setbacks = block_data.push_features(special_setbacks, features = 'B01001e1')

    print('Special lots population', special_lots['B01001e1'].sum())
    print('Special setbacks population', special_setbacks['B01001e1'].sum())

def analyze_landmarks_in_houston():

    # Get points and polygons
    houston_landmarks = gpd.read_file(zoning.houston_historic_landmarks_path).to_crs({'init': 'epsg:4326'})
    tx_hd_data = gpd.read_file(zoning.tx_hd_path)
    houston_nat_districts = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston'].to_crs({'init': 'epsg:4326'})
    block_data = boundaries.BlockBoundaries(['X01_AGE_AND_SEX'], cities='Houston', get_percent_residential = False)

    # Calculate number of points in polygons
    num_points = utilities.spatial_joins.points_intersect_multiple_polygons(houston_landmarks, houston_nat_districts).count()
    houston_nat_districts = block_data.push_features(houston_nat_districts, features = 'B01001e1', account_method='water')
    population = houston_nat_districts['B01001e1'].sum()

    print('Number of Historic Landmarks in NRHDS is {}, Population of NRHDs is {}'.format(num_points, population))

def analyze_austin_overlays():
    """Measures area in Austin conditional on distance from city center and number of zoning overlays"""

    austin_zones = zoning.austin_inputs.process_zoning_shapefile()

    # Isolate the overlay by excising the base zone
    def get_overlay(row):
        row['overlay'] = str(row['zoning_zty']).replace(str(row['base_zone']), '')
        return row
    austin_zones = austin_zones.apply(get_overlay, axis = 1)
    austin_zones.crs = {'init':'epsg:4326'}

    # Get other important features, i.e. area, dist to center
    austin_zones = utilities.measurements.get_area_in_units(austin_zones)
    austin_zones['dist_to_center'] = utilities.measurements.calculate_dist_to_center(austin_zones,
                                                                                     zoning.austin_inputs.lat,
                                                                                     zoning.austin_inputs.long)

    austin_zones['dist_to_center'] = austin_zones['dist_to_center'].apply(lambda x: 2*np.ceil(x/2))

    # Calculate the number of overlays (each overlay is separated by a '-')
    austin_zones['num_overlays'] = austin_zones['overlay'].apply(lambda x: str(x).count('-'))

    # Find areas conditional on dist_to_center and number of overlays, then normalize to add up to 100% for each
    # value for dist_to_center
    result = austin_zones.groupby(['dist_to_center', 'num_overlays'])['area'].count()
    divisor = austin_zones.groupby(['dist_to_center'])['area'].count()
    result = result.divide(divisor)
    print('Here are the areas in Austin conditional on distance to city center and number of zoning overlays:')
    print(result)