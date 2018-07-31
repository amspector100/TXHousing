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

def process_harris_parcel_data():

    time0 = time.time()
    print('Reading Houston parcel data')
    parcel_data = gpd.read_file(harris_parcel_path_2018)
    print('Finished reading Houston parcel data, time is {}. Now transforming.'.format(time.time() - time0))

    # Get rid of bad geometries and process crs. This eliminiates 5,000 out of nearly 1.5 million parcels so is rather negligible.
    parcel_data = parcel_data.loc[~parcel_data['geometry'].apply(lambda x: x is None)]
    parcel_data = parcel_data.loc[~parcel_data['geometry'].apply(lambda x: x.is_empty)]
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


    return houston_parcels
