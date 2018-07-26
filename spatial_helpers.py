import numpy as np
import datetime as dt
import time
import pandas as pd
import geopandas as gpd
import fiona
import shapely.geometry
import matplotlib.pyplot as plt
from tqdm import tqdm

from inputs import *
import spatial_functions as sf

# This is a bit of a weird file. It's for helper (data processing) functions that require spatial transformations
# and therefore require importing the spatial_functions model. Most of these functions are run once and then cached and
# never run again (hopefully).

def process_travis_parcel_data():

    time0 = time.time()
    print('Reading Houston parcel data')
    parcel_data = gpd.read_file(travis_parcel_path_2018)
    print('Finished reading Houston parcel data, time is {}. Now transforming.'.format(time.time() - time0))

    # Get rid of bad geometries and process crs. This eliminiates 5,000 out of nearly 1.5 million parcels so is rather negligible.
    parcel_data = parcel_data.loc[[not bool for bool in parcel_data['geometry'].apply(lambda x: x is None)]]
    parcel_data = parcel_data.loc[[not bool for bool in parcel_data['geometry'].apply(lambda x: x.is_empty)]]
    parcel_data = parcel_data.loc[parcel_data['geometry'].apply(lambda x: isinstance(x, shapely.geometry.polygon.Polygon))]
    parcel_data = parcel_data.to_crs({'init':'epsg:4326'})
    print('Finished transforming, time is {}. Now subsetting.'.format(time.time() - time0))

    # Reset index
    parcel_data.reset_index(inplace = True)

    # Get houston shape
    place_shapes = gpd.read_file(texas_places_path)
    houston_shape = place_shapes.loc[place_shapes['NAME'] == 'Houston']
    houston_grid = sf.fragment(houston_shape['geometry'].iloc[0], 12, 12) # List of polygons - not a geodataframe

    # Intersection by centroids for efficiency
    parcel_data['centroids'] = parcel_data['geometry'].centroid
    centroids = parcel_data[['centroids']].set_geometry('centroids')

    # Use grid and sindex to find intersections
    spatial_index = centroids.sindex
    true_intersections = set()
    for grid_piece in tqdm(houston_grid):
        possible_intersections_index = list(spatial_index.intersection(grid_piece.bounds))
        possible_intersections = centroids.iloc[possible_intersections_index]
        precise_intersections_bools = possible_intersections['centroids'].intersects(grid_piece)
        precise_intersections = set(possible_intersections[precise_intersections_bools].index)
        true_intersections = true_intersections.union(precise_intersections)
    true_intersections = list(true_intersections)
    print('Finished subsetting. Time is {}. Now saving and plotting for accuracy.'.format(time.time() - time0))

    # Only consider ones in houston
    houston_parcels = parcel_data.loc[true_intersections]
    houston_parcels = houston_parcels[[col for col in houston_parcels.columns if col != 'centroids']]

    print('Saving')
    houston_parcels.to_file(houston_parcel_path_2018)
    print('Plotting')
    houston_parcels.plot()
    plt.show()

    return houston_parcels

# Hist of property values of downzoning and construction vs the city at large for Austin specifically.
def process_austin_parcel_data(zip_intersections = True, save = True):
    """
    Processes austin parcel data. Normalizes property values by land use, computes zip code intersections if zip_intersections = True,  and formats base
    zone properly, and then returns. Will save as a shapefile if save = True.
    :return:
    """

    time0 = time.time()
    print('Starting to read parcel data')
    parcel_data = gpd.read_file(inputs.austin_parcel_path)
    print('Finished reading parcel data, took {} seconds. Now processing base zones.'.format(time.time() - time0))

    # Get basezone in proper format
    def format_basezone(s):
        s = str(s)
        try:
            return s.split('(')[1][:-1]
        except:
            # Usually this means there's no basezone information
            return s
    parcel_data['basezone'] = parcel_data['basezone'].apply(format_basezone)

    # Compute value per area and average value per area for each land use, so we can normalize later
    print('Finished processing base zones, took {} seconds. Now normalizing data conditionally on land use.'.format(time.time() - time0))
    parcel_data.loc[:, 'value_per_area'] = parcel_data['appraised'].divide(parcel_data['shape_area'])


    land_use_means = {}
    land_use_stds = {}
    for code in parcel_data['land_use'].unique():
        subset = parcel_data.loc[parcel_data['land_use'] == code]
        land_use_means[code] = subset['value_per_area'].mean()
        land_use_stds[code] = subset['value_per_area'].std()

    # Normalize - there should be no key errors because we just called parcel_data['land_use'].unique().
    vpa_mean = parcel_data['land_use'].map(land_use_means)
    vpa_std = parcel_data['land_use'].map(land_use_stds)
    norm_pva = (parcel_data['value_per_area'] - vpa_mean).divide(vpa_std)
    parcel_data['normalized_value_per_area'] = norm_pva


    # Zip intersections
    if zip_intersections:
        print('Finished normalizing, took {} seconds. Computing zip code intersections.'.format(
            time.time() - time0))
        parcel_data = sf.zip_intersect(parcel_data, zips = inputs.austin_zips)


    if save:
        print('Finished processing parcel data, took {} seconds. Now writing to shapefile.'.format(time.time() - time0))
        parcel_data.to_file(austin_processed_parcel_data_path)
        print('Finished writing to shapefile, took {} seconds.'.format(time.time() - time0))
    else:
        print('Finished processing parcel data, took {} seconds.'.format(time.time() - time0))

    return parcel_data
