import numpy as np
import datetime as dt
import time
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry.polygon import *
from shapely.geometry.multipolygon import *
from shapely.geometry.multipoint import *
from shapely.geometry.point import *
import shapely
from shapely.ops import cascaded_union, polygonize, unary_union
from shapely import geometry
from scipy.spatial import Delaunay
import math
from tqdm import tqdm
import helpers
import copy
from math import radians, cos, sin, asin, sqrt, pi
import pyproj
from inputs import *
from collections import Iterable


# Everything in this file should be sufficiently general as to work without modifications for any city in the United States

# Utilities ------------------------------------------------------------------------------------------------------------

# Get urban cores
def get_urban_core(zoning_input, radius, scale = 5280, newproj = 'epsg:2277'):
    """
    Create a polygon representing the urban core of a city. Basically this is the buffer function of shapely but it uses crs transformations so you get to pick the units.
    :param zoning_input: A zoning input, used to find the center of the city
    :param radius: The radius in miles or km
    :param scale: Defaults to 5280, feet per mile.
    :param newproj: The new projection to use to calculate this distance (by default epsg:2277, which is in feat).
    :return: a geopandas geodataframe with a single column (geometry) of length one (polygon) which represents the urban core.
    """
    # Create point, transform to new coords. Do long lat because there's no standardization.
    core = gpd.GeoDataFrame(geometry = [shapely.geometry.point.Point(zoning_input.long, zoning_input.lat)])
    core.crs = {'init':'epsg:4326'}
    core = core.to_crs({'init':newproj})

    # Get point and create buffer
    core = core['geometry'][0]
    core = core.buffer(scale*radius)
    core = gpd.GeoDataFrame(geometry = [core])
    core.crs = {'init':newproj}

    # Return to dataframe and transform back to lat long
    core = core.to_crs({'init':'epsg:4326'})
    return core

# Get area in miles. GDF must initially be in lat long
def get_area_in_units(gdf, geometry_column = 'geometry', newproj = 'epsg:2277', scale = 3.58701*10**(-8), name = 'area'):
    """
    Get area of polygons of a geodataframe in units. By default, gets it in miles.
    :param gdf: Geodataframe with polygons in the geometry column.
    :param geometry_column: Geometry column of the geodataframe, defaults to 'geometry'
    :param newproj: The new projection to use to calculate units. Defaults to epsg:2277, which is probably fine for
    Austin/Dallas/Houston and is in feet.
    :param scale: A scale to multiply by. Defaults to 3.58701*10**(-8) which is the number of square miles in a square foot.
    :param name: the name of the new column that will be created to store the area information. Defaults to 'area'.
    :return: The geodataframe with a column named name (defualts to 'area') which has the area of each polygon in the desired units.
    """
    old_projection = gdf.crs
    gdf = gdf.to_crs({'init':newproj})
    gdf['area'] = scale*gdf[geometry_column].area
    gdf = gdf.to_crs(old_projection)
    return gdf

# Block data ------------------------------------------------------------------------------------------------------------------block data

# See datadic at https://www2.census.gov/geo/tiger/TIGER_DP/2016ACS/Metadata/BG_METADATA_2016.txt
def get_block_geodata(data_layers, cities=None, get_percent_residential = True):
    """
    Get geodata by block group and subset to only include the municipality
    :param data_layers: Iterable of codes for the data layer of the geodatabase.
    :param cities: City name (i.e. 'Austin'), or iterable of city names to filter by (i.e. ['Austin', 'Dallas'].
    Need to be in Texas or this won't work. Defaults to None, and will not filter the data at all.
    :return: Dictionary of geodataframes in the form {cityname: geodataframe} or just a single geodataframe if 'cities'
    is a string, not an iterable.
    """

    time0 = time.time()
    print('Getting block geodata for {}'.format(cities))
    # Merge layers from geodatabase - this is a pretty cheap operation, relatively speaking (10 sec or so)
    geodata = gpd.read_file(texas_blocks_path, layer='ACS_2016_5YR_BG_48_TEXAS')
    geodata['geometry'] = geodata['geometry'].apply(lambda x: x[0])
    for data_layer in data_layers:
        data = gpd.read_file(texas_blocks_path, layer=data_layer)
        geodata = geodata.merge(data, how='inner', left_on='GEOID_Data', right_on='GEOID')
        geodata.rename({'GEOID_y': 'GEOID'}, axis='columns', inplace=True)

    geodata.set_index('GEOID', inplace = True)

    # Get the percent of land which is zoned residential inside the city limits
    if get_percent_residential:
        percent_residential = pd.read_csv('data/bg_percent_residential.csv', index_col = 0).fillna(1)
        percent_residential.columns = ['percent_residential']
        geodata = geodata.join(percent_residential)

    geodata = gpd.GeoDataFrame(data = geodata[[x for x in geodata.columns if x != 'geometry_x']], geometry = geodata['geometry_x'])
    geodata.rename({'geometry_x': 'geometry'}, axis = 'columns', inplace = True)
    geodata.crs = {'init':'epsg:4326'}
    spatial_index = geodata.sindex

    if cities is not None:

        # Only consider blocks in cities
        place_shapes = gpd.read_file(texas_places_path)

        def fast_intersect(city):
            shape = place_shapes.loc[place_shapes['NAME'] == city, 'geometry'].values[0]
            possible_intersections_index = list(spatial_index.intersection(shape.bounds))
            possible_intersections = geodata.iloc[possible_intersections_index]
            precise_intersections = possible_intersections['geometry'].intersects(shape)
            return precise_intersections[precise_intersections.values] # only include shapes that actually intersect

        # Check if cities is a string, if so get a result
        if isinstance(cities, str):
            result = geodata.loc[fast_intersect(cities).index.tolist()]
        else:
            result = {}
            for city in cities:
                precise_intersections = fast_intersect(city)
                result[city] = geodata.loc[precise_intersections.index.tolist()]
        print('Finished getting block geodata, time is {}'.format(time.time() - time0))
        return result

    else:
        print('Finished getting block geodata, time is {}'.format(time.time() - time0))
        return geodata

def get_average_by_area(data_source, spatial_index, polygon, density_features = 'B01001e1', geometry_column = 'geometry'):
    """
    Calculates the average 'feature' of a 'polygon' using a 'data_source' of different shapes which (in some combination)
    cover the polygon. If you want to do this for a list of polygons, use get_all_averages_by_area, which is listed below.
    :param data_source: The data source, usually block data. Geodataframe. Must have a column specifying total population.
    :param spatial_index: Spatial index of data_source. You can generate this by calling 'data_source.sindex' right
    before calling this function.
    :param polygon: A polygon to find the population density of. Note that this should be a polygon, not a geopandas
    object.
    :param density_features: List of features (or just a single feature as a string) in the data_source of interest.
     Defaults to 'B01001e1' which is the total population est in block data 'X01_AGE_AND_SEX' layer. Important: This
    feature MUST be of the form 'units/area' where the area units (i.e. square feet) are the same as the units of the
    area in the geopandas .area attribute.
    :param geometry_column: Geometry column of the block data, defaults to 'geometry'
    :return: float (the feature of the area, weighted by area)
    """

    # If the feature is a string, turn it to a list
    if isinstance(density_features, str):
        density_features = [density_features]

    # Find which blocks intersect
    possible_intersections_index = list(spatial_index.intersection(polygon.bounds))
    possible_intersections = data_source.iloc[possible_intersections_index]
    precise_intersections_index = possible_intersections[geometry_column].intersects(polygon)
    precise_intersections = data_source.loc[precise_intersections_index[precise_intersections_index.values].index.tolist()]

    # Now actually find the intersections. Possible buffer them to avoid absurdly high numbers (area can be super small).
    precise_intersections.loc[:, geometry_column] = precise_intersections[geometry_column].intersection(polygon)
    precise_intersections.loc[:, 'area'] = precise_intersections[geometry_column].area
    precise_intersections = precise_intersections.loc[precise_intersections['area'] != 0] # For some reason we get 0 area every now and then, so get rid of these columns.

    # Now find feature for the polygon
    result = precise_intersections[density_features].transpose().dot(precise_intersections['area'])
    return result

def get_all_averages_by_area(data_source, other_geometries, features = 'B01001e1', data_source_geometry_column = 'geometry',
                             other_geometries_column = 'geometry', drop_multipolygons = True, fillna = None, account_method = 'percent_residential'):
    """
    Get averages of a 'feature' from 'data_source' by area. Note: data_source and other_geometries should have the
    same crs initially.
    :param data_source: The data source, usually block data. Geodataframe. Must have a column specifying total population.
    :param other_geometries: Geodataframe, where the geometry column is filled with polygons. Will calculate the feature
    for each of these polygons.
    :param features: The feature in question, defaults to 'B01001e1' which is the population estimate in the X01_SEX_AND_AGE
    layer in block geodata. Can also be a list of features, i.e. ['B01001e1', 'B01001e2']
    :param density: Boolean. Default False. If True, will assume that the 'feature' is in units per area and will not
    divide the feature by the area of the data source polygons.
    :param data_source_geometry_column: geometry column for data_source
    :param other_geometries_column: geometry column for other_geometries
    :param fillna: if not None, fill na values with this value.
    :param account_method: The method by which to account for the % of an area which is not residential (this prevents
    population-related estimates from being too low). Can either be 'None', 'percent_residential', or 'percent_land'
    :return: other_geometries but with a new column, feature, which has the averages by area.
    """
    time0 = time.time()
    if isinstance(features, str):
        features = [features]

    if data_source.crs != other_geometries.crs:
        print("""Note: in get_all_averages_by_area, the crs for data_source ({})  and other_geometries ({}) is not the same. This
              'might cause issues.""".format(data_source.crs, other_geometries.crs))

    # Process data for convenience (just to prevent multipolygons/invalid polygons from messing things up)
    data_source = process_geometry(data_source, drop_multipolygons=drop_multipolygons)
    data_source.reset_index(drop = True)
    data_source.index = [str(ind) for ind in data_source.index]

    # Get the feature in terms of units per area (do a bit of renaming to make it clear these are densities)
    old_columns_dictionary = {str(feature) + '_density':feature for feature in features}
    new_columns_dictionary = {feature:str(feature) + '_density' for feature in features}
    new_columns = [new_columns_dictionary[key] for key in new_columns_dictionary]
    if account_method == 'water':
        data_source.loc[:, 'percent_land'] = data_source['ALAND'].divide(data_source["ALAND"] + data_source['AWATER'])
        densities = data_source[features].multiply(data_source['percent_land'], axis = 0).divide(data_source[data_source_geometry_column].area, axis = 0)
    elif account_method == 'percent_residential':
        data_source = data_source.loc[data_source['percent_residential'] != 0]
        densities = data_source[features].divide(data_source['percent_residential'], axis = 0).divide(data_source[data_source_geometry_column].area, axis = 0)
    else:
        densities = data_source[features].divide(data_source[data_source_geometry_column].area)

    # Rename and join to data
    densities = densities.rename(columns = new_columns_dictionary)
    data_source = data_source.join(densities)

    other_geometries = process_geometry(other_geometries, drop_multipolygons = drop_multipolygons)
    if len(other_geometries) == 0 and drop_multipolygons == True:
        # Warn user in case
        print('Warning in get_all_averages_by_area: it looks like dropping multipolygons in the'
              '"process geometry" call has eliminated all the data.')

    # Get spatial index
    spatial_index = data_source.sindex

    # Quick function to apply to geometry column for other_geometries
    def get_avg(polygon):
        result = get_average_by_area(data_source, spatial_index, polygon, density_features = new_columns, geometry_column = data_source_geometry_column)
        if fillna is not None:
            result.fillna(fillna, inplace = True)
        return result

    # Get averages by area - this takes a while.
    final_values = other_geometries[other_geometries_column].apply(get_avg)
    final_values = final_values.rename(columns = old_columns_dictionary)
    other_geometries = other_geometries.join(final_values)
    print("Finished calculating average {} by area, took {}".format(features, time.time() - time0))

    return other_geometries





# Fyi, list of block layers from geodatabase is as follows:
# ['X00_COUNTS', 'X01_AGE_AND_SEX', 'X02_RACE', 'X03_HISPANIC_OR_LATINO_ORIGIN', 'X07_MIGRATION', 'X08_COMMUTING',
# 'X09_CHILDREN_HOUSEHOLD_RELATIONSHIP', 'X11_HOUSEHOLD_FAMILY_SUBFAMILIES', 'X12_MARITAL_STATUS_AND_HISTORY',
# 'X14_SCHOOL_ENROLLMENT', 'X15_EDUCATIONAL_ATTAINMENT', 'X16_LANGUAGE_SPOKEN_AT_HOME', 'X17_POVERTY', 'X19_INCOME',
# 'X20_EARNINGS', 'X21_VETERAN_STATUS', 'X22_FOOD_STAMPS', 'X23_EMPLOYMENT_STATUS', 'X24_INDUSTRY_OCCUPATION',
# 'X25_HOUSING_CHARACTERISTICS', 'X27_HEALTH_INSURANCE', 'X99_IMPUTATION', 'BG_METADATA_2016', 'ACS_2016_5YR_BG_48_TEXAS']

# In X01_AGE_AND_SEX, the following are useful:
# B01001e1: SEX BY AGE: Total: Total population -- (Estimate)


# In X19_INCOME, the following are useful:
# B19013e1:	MEDIAN HOUSEHOLD INCOME IN THE PAST 12 MONTHS (IN 2016 INFLATION-ADJUSTED DOLLARS) (Estimate)


# Basically meant to simplify very complex zoning data for categorical features ---------------------------------------------------------------------------------------------------------

# Global granularity variable
max_comb = 30
alpha = 0.01 # I got this by experimenting

# Check if something is a multipolygon
def is_polygon(something):
    try:
        result = something.type == 'Polygon'
        return result
    except:
        return False

# Make valid by buffering
def make_valid_buffer(poly):
    if poly.is_valid == False:
        return poly.buffer(0)
    else:
        return poly

# Ignore multipolygons and make invalid polygons valid
def process_geometry(gdf, geometry_column = 'geometry', drop_multipolygons = True):
    gdf.loc[:, geometry_column] = gdf[geometry_column].apply(make_valid_buffer)
    if drop_multipolygons:
        gdf = gdf.loc[gdf[geometry_column].apply(is_polygon)]
    return gdf

# Drop invalid points
def process_points(points, geometry_column = 'geometry'):
    points = points.loc[points[geometry_column].is_valid]  # Ignore invalid points (i.e. with 'na's)
    points.reset_index(drop=True, inplace=True)
    points.index = [str(ind) for ind in points.index]
    return points


# Adapted from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
def alpha_shape(points, alpha = 0.001):
    """
    Compute the alpha shape/concave hull of points.
    :param points: Iterable container of shapely points.
    :param alpha: alpha value to influence the the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """

    # This algorithm is pointless unless there are more than 3 points
    if len(points) < 4:
        return MultiPoint(list(points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between points i and j if not already in the list
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    # get coordinates
    coords = np.array([point.coords[0]for point in points])

    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # Loop through triangles
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0
        # Area of triangle by Heron's formula. Sometimes this comes out as a negative number, who knows why.
        try:
            area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        except:
            continue
        # Prevent div by 0 errors
        if area == 0:
            continue
        # Get circum
        circum_r = a * b * c / (4.0 * area)
        # Here's the radius filter.
        # print circum_r
        if circum_r < 1.0 / alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    result = cascaded_union(triangles)
    if result.is_valid == False:
        result = result.buffer(0)
    return result # could return edge_points

# Convert multipolygon to polygon - this doesn't have to be fast, will be used on at most like 0.5% of the objects
def multipolygon_to_polygon(mp):
    parts = list(mp)
    try:
        union = cascaded_union(parts)
        return union
    except:
        points = []
        for part in parts:
            points.extend([Point(coord) for coord in part.exterior.coords])
        return alpha_shape(points, alpha = alpha)

# Combine data
def combine(data, feature, spatial_index, k, max_comb = max_comb, alpha = alpha):

    # Get polygon and zone
    test_polygon = data.loc[k, 'geometry']
    test_zone = data.loc[k, feature]

    # Create index and find polygons
    new_index = spatial_index.nearest(test_polygon.bounds, num_results = max_comb) # .16 milliseconds on average
    polys = data.loc[list(new_index)]
    good_polys = [Point(coord) for coord in list(test_polygon.exterior.coords)]
    # good_polys =
    good_inds = [k]

    # Loop through to see which ones are in the right zone
    for ind, row in polys.iterrows():
        if row[feature] == test_zone and row['flag'] == False:
            good_polys.extend([Point(coord) for coord in list(row['geometry'].exterior.coords)])
            #good_polys.append(row['geometry'])
            good_inds.append(ind)
        else:
            break

    # Get convex hull and return
    #result = MultiPolygon(good_polys).convex_hull
    result = alpha_shape(good_polys, alpha)# This is the part that takes a while
    return result, good_inds, test_zone

# Combine data - new idea - unfinished
def combine_v2(data, feature, spatial_index, k, max_comb = max_comb, alpha = alpha):

    # Get polygon and zone
    test_polygon = data.loc[k, 'geometry']
    test_zone = data.loc[k, feature]

    # Create index and find polygons
    new_index = spatial_index.nearest(test_polygon.bounds, num_results = max_comb) # .16 milliseconds on average
    polys = data.loc[list(new_index)]

    # Possible matches
    possible_polys = polys.loc[polys['broad_zone'] == test_zone]
    # Non-matches
    bad_polys = polys.loc[polys['broad_zone'] != test_zone]
    # Initialize list of good matches, both indices and coordinates
    good_poly_coords = []
    good_inds = [k]

    for ind, row in possible_polys.iterrows():
        possible_poly = row['geometry']
        line = geometry.linestring.LineString([test_polygon.centroid.coords[:][0], possible_poly.centroid.coords[:][0]])
        if any([line.intersects(bad_poly) for bad_poly in bad_polys['geometry']]):
            continue
        else:
            good_poly_coords.extend([Point(coord) for coord in list(row['geometry'].exterior.coords)])
            good_inds.append(ind)

    result = alpha_shape(good_poly_coords, alpha)# This is the part that takes a while
    return result, good_inds, test_zone

def combine_all(data, feature, max_comb = max_comb, centroids = False, ignore_features = False, to_ignore = ['Other'],
                alpha = alpha, use_v2 = False):
    """
    :param data: Geodataframe. Index should be a rangeindex. Geometry should be polygons. One other column.
    :param max_comb: Granularity.
    :param use_v2: Use a slightly different combine tactic.
    :return: Geodataframe with a rangeindex and polygon geometry (as well as that one other column).
    """

    # Start timing
    time0 = time.time()

    # Spatial index, either by centroids or by polygons
    if centroids:
        data['centroids'] = data['geometry'].centroid
        data = data.set_geometry('centroids')
    spatial_index = data.sindex

    # If the flag is 0, then the polygon in that row has not been combined yet.
    all_good_inds = []
    data['flag'] = 0
    result = gpd.GeoDataFrame(pd.DataFrame(columns = [feature, 'geometry']), crs = {'init':'EPSG:4326'})

    # Run through
    for i in tqdm(np.arange(0, len(data), 1)):

        # If it's already been dealt with, ignore it
        try:
            if ignore_features == True and data.loc[i, feature] in to_ignore:
                continue
            elif data.loc[i, 'flag'] == 1:
                continue
        except KeyError as e:
            print('Almost got {} error'.format(e))
            print('Here is the index: {}'.format(data.index))
            print('i = {}'.format(i))

        # Combine polygons around it
        if use_v2:
            poly, good_inds, zone = combine_v2(data, feature, spatial_index, i, max_comb = max_comb, alpha = alpha)
        else:
            poly, good_inds, zone = combine(data, feature, spatial_index, i, max_comb = max_comb, alpha = alpha)

        data.at[good_inds, 'flag'] = 1
        all_good_inds.extend(good_inds)
        result.loc['new' + str(i)] = [zone, poly]

    # Reindex for convenience
    result.index = np.arange(0, len(result.index), 1)

    print('Data of length {} took {}, cut it down to length {}'.format(len(data), time.time() - time0, len(result)))
    return(result)

def process_combined_result(data, factor = 'broad_zone'):
    """
    Takes unary unions of each type of zone or feature to reduce the size of the result. Also makes sure nothing overlaps.
    :param result: Result of a combine_all call
    :return: gdf with one row per unique value in the factor column.
    """

    result = gpd.GeoDataFrame(columns = ['geometry'])
    for zone in data[factor].unique():
        result.at[zone, 'geometry'] = data.loc[data[factor] == zone, 'geometry'].unary_union
    result['area'] = result['geometry'].area
    result.sort_values(by = 'area', ascending = True, inplace = True)

    print('Combining final result of combine_all call')
    for zone in tqdm(result.index):
        if zone == result.index[0]:
            running_union = result.loc[zone, 'geometry']
        else:
            cached_geometry = copy.copy(result.loc[zone, 'geometry'])
            result.at[zone, 'geometry'] = unary_union(cached_geometry.difference(running_union))
            running_union = cached_geometry.union(running_union)

    result[factor] = result.index
    return result

# Intersect zoning with zip codes --------------------------------------------------------------------------------------------------------------------------------------
def zip_intersect(gdf, zips):
    """
    :param gdf: A geodataframe of some sort (it should probably be in the US, otherwise this is pointless).
    :param zips: The zip codes in the city.
    :return: The gdf but with an extra column for each zip code which specifies the fractional area in the zip code.
    """

    # Make sure zips input and data index is full of strings
    zips = [str(something) for something in zips]
    gdf.index = [str(ind) for ind in gdf.index]
    gdf = gdf.to_crs({'init': 'EPSG:4326'})

    # Get the zip code geodata
    zipdata = helpers.get_zip_boundaries()
    zipdata = zipdata.loc[zipdata['ZCTA5CE10'].isin(zips), ['geometry']]

    # Create spatial index of zoning data
    spatial_index = gdf.sindex

    # Run through zip codes
    print('Starting to calculate which zones are in which zip codes')
    for code, row in tqdm(zipdata.iterrows()):

        # Find intersections
        polygon = row['geometry']
        possible_matches_index = list(spatial_index.intersection(polygon.bounds))
        possible_matches = gdf.iloc[possible_matches_index]
        precise_matches = possible_matches['geometry'].intersection(polygon)

        # Calculate percent of area in the zip code - this is a bit wasteful because most zones are only in one zip code,
        # but it avoids fragmenting the shapes which are quite hard to put back together.
        gdf[code] = precise_matches.area / gdf['geometry'].area
        gdf[code].fillna(value = 0, inplace = True)

    return gdf


# Given a bunch of points, it splits them into a rectangular grid of m by n and calculates either the number of points
# in each rectangle or the mean or median of a factor associated with the points for each rectangle. ---------------------------------------------------------------------

# Helper function for spatial tree efficiency - fragment polygon into smaller sizes
def fragment(polygon, horiz, vert):
    """
    :param polygon: Polygon to fragment
    :param horiz: # of horizontal fragments
    :param vert: # of vertical fragments
    :return: Grid, a list of polygons
    """
    minx, miny, maxx, maxy = polygon.bounds
    xlen = (maxx - minx) / horiz
    ylen = (maxy - miny) / vert
    grid = []
    for i in np.arange(0, horiz, 1):
        for j in np.arange(0, vert, 1):
            b = shapely.geometry.box(xlen * i + minx, ylen * j + miny, xlen * (i + 1) + minx,
                    ylen * (j + 1) + miny)  # Left, botton, right, upper
            g = polygon.intersection(b)
            if g.is_empty:
                continue
            grid.append(g)

    return grid

def make_point_grid(gdf, horiz=20, vert=20, factor=None, by='mean', geometry_column='geometry'):
    """
    :param gdf: Geodataframe, with point geometry presumably.
    :param horiz: Horizontal number of boxes
    :param vert: Vertical number of boxes
    :param factor: Defaults to None, otherwise the (continuous) value with which to take the mean/median of the points
    for the choropleth color.
    :param by: 'mean' or 'median'. Meaningless unless you have the factor column.
    :param geometry_column: The column the points are contained in. These should be shapely points in lat/long form.
    :return: geodataframe with grid geometry and a 'value' column
    """

    # Work with crs's for a sec
    if gdf.crs is None:
        gdf.crs = {'init': 'epsg:4326'}

    # Adapted from snorfalorpagus - could make this a helper function in the future if necessary
    points = shapely.geometry.multipoint.MultiPoint(gdf[geometry_column].tolist())
    hull = points.convex_hull
    #minx, miny, maxx, maxy = points.bounds
    #xlen = (maxx - minx) / horiz
    #ylen = (maxy - miny) / vert
    #grid = []
    #for i in np.arange(0, horiz, 1):
        #for j in np.arange(0, vert, 1):
            #b = shapely.geometry.box(xlen * i + minx, ylen * j + miny, xlen * (i + 1) + minx,
            #        ylen * (j + 1) + miny)  # Left, botton, right, upper
            #g = hull.intersection(b)
            #if g.is_empty:
            #    continue
            #grid.append(g)

    grid = fragment(hull, horiz = horiz, vert = vert)

    # Now initialize result
    result = gpd.GeoSeries(grid)
    result.index = [str(ind) for ind in result.index]  # prevent ugliness and silent errors arising from spatial index
    value = pd.Series(index=result.index)

    # Get intersections and fill values
    spatial_index = gdf.sindex
    for ind in result.index:
        box = result[ind]
        possible_matches_index = list(spatial_index.intersection(box.bounds))
        possible_matches = gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[geometry_column].intersection(box)

        if factor is None:
            value[ind] = len(precise_matches)
        elif by == 'mean':
            value[ind] = precise_matches.loc[precise_matches[factor].notnull(), factor].mean()
        elif by == 'median':
            value[ind] = precise_matches.loc[precise_matches[factor].notnull(), factor].median()
        else:
            print(
            'In point_choropleth call, "by" argument must either equal "mean" or "median" - you put "{}"'.format(by))

    result = gpd.GeoDataFrame(data = pd.DataFrame(data = value.values, columns = ['value'], index = value.index), geometry = result)
    result.crs = gdf.crs
    return result

# Create dist from city center graphs --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Calculates distance between two lat/long points using haversine formula.
# Reference https://gis.stackexchange.com/questions/279109/calculate-distance-between-a-coordinate-and-a-county-in-geopandas
def haversine(point1, point2):
    """
    :param point1: Shapely point. Long then lat.
    :param point2: Shapely point. Long then lat.
    :return: Distance in miles.
    """

    lon1, lat1 = list(point1.coords[:][0][0:])
    lon2, lat2 = list(point2.coords[:][0][0:])

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 3956  # Radius of earth in miles. Use 6371 for km
    return c * r


def points_intersect_rings(gdf, zoning_input, factor = None, step = 1, categorical = True, by = 'mean',
                           geometry_column = 'geometry', per_square_mile = True):

    """
    :param gdf: Geodataframe. Assumes this is already transformed to lat long coord system. Should be polygon geometry.
    :param factor: A factor to calculate for. If none, will just calculate the area or number of points radiating out from city center.
    :param zoning_input: Zoning input, at the center of which we start our rings
    :param step: Number of miles where the ring radiates outwards.
    :param maximum: Max radius (miles)
    :param categorical: If true, will only calculate percent land (or % of points) in the radius
    :param by: Defaults to "mean". If categorical = False, use "by" to determine how to calculate averages over points.
    :param geometry_column: name of geometry column. Default geometry.
    :param psm: if true, divide by the area of the ring.
    :return: Dataframe
    """

    # Get distance from center of city, in miles currently
    center = shapely.geometry.point.Point(zoning_input.long, zoning_input.lat)

    def rounded_dist_to_center(point):
        dist = haversine(point, center)
        rdist = step*((dist // step) + 1) # Do this to get smoother bins (i.e. always round up)
        return rdist

    gdf['dist_to_center'] = gdf[geometry_column].apply(rounded_dist_to_center)

    # Get counts in this case
    if factor is None:
        result = gdf[['geometry', 'dist_to_center']].groupby('dist_to_center')['geometry'].count()

    # Get counts by factor level
    elif categorical:
        result = gdf[['dist_to_center', factor]].groupby(['dist_to_center', factor]).size()
        result = result.unstack().fillna(0)

    # Get mean and median
    elif by == 'mean':
        result = gdf[[factor, 'dist_to_center']].groupby('dist_to_center')[factor].mean()
    elif by == 'median':
        result = gdf[[factor, 'dist_to_center']].groupby('dist_to_center')[factor].median()

    # Warn user if they supplied a bad arg
    else:
        print("""In points_intersect_rings call, when categorical = False and factors is not None, "by" argument must
        either equal "mean" or "median" - you put "{}" """.format(by))
        result = None

    # If told, get area per square mile
    if per_square_mile:

        def get_area_ring(radius):
            return pi*(radius**2 - (radius - step)**2)

        areas = pd.Series(result.index.map(get_area_ring).tolist(), index = result.index)

        if isinstance(result, pd.DataFrame):
            result = result.divide(areas, axis = 0)
        else:
            result = result.divide(areas)

    # Return result
    return result


# Radius from city center for polygon data
def polygons_intersect_rings(gdf, zoning_input, factor = None, newproj = 'epsg:2277', step = 1, maximum = 20,
                    categorical = True, geometry_column = 'geometry', city = None):

    """
    :param gdf: Geodataframe. Assumes this is already transformed to lat long coord system. Should be polygon geometry.
    :param factor: A factor to calculate for. If none, will just calculate the area or number of points radiating out from city center.
    :param zoning_input: Zoning input, at the center of which we start our rings
    :param newproj: the new projection system necessarily used in this. Defaults to 2277 which is good for Austin and fine for Texas. Note units in this are in feet.
    :param step: Number of miles where the ring radiates outwards.
    :param maximum: Max radius (miles)
    :param points: Boolean. If true, assumes geometry of gdf is a point geometry.
    :param geometry_column: name of geometry column. Default geometry.
    :param categorical: If true, will only calculate percent land (or % of points) in the radius. Else will calculate mean by area.
    :param city: If not "none", if factor is "none", will read the shaepfile of the boundaries of the city to ensure
    more accurate calculations. (Otherwise, for a radius of size 12, the area of the circle might be greater than the area
    of the city).
    :return: Dataframe or Series
    """

    feet_to_mile = 5280
    gdf = process_geometry(gdf)
    gdf.reset_index(drop = True)
    gdf.index = [str(ind) for ind in gdf.index]


    # Necessarily will have to transform into new epsg code in order to efficiently calculate this.
    # Here, gdf should be in (lat, long) format
    if gdf.crs is None:
        print('No CRS set, assuming lat and long data before transformation.')
        gdf.crs = {'init':'epsg:4326'}

    gdf = gdf.to_crs({'init':newproj})


    # Get center and initialize output
    center_latlong = shapely.geometry.point.Point(zoning_input.long, zoning_input.lat)
    center_series = gpd.GeoSeries(center_latlong)
    center_series.crs = {'init':'epsg:4326'}
    center_gdf = gpd.GeoDataFrame(geometry = center_series)
    center_gdf = center_gdf.to_crs({'init':newproj})
    center = center_gdf.loc[0, 'geometry']


    # Initialize result. If categorical, need a dataframe (one column for each unique value). Else, use a pd.Series.
    if factor is None or categorical == False:
        result = pd.Series()
        if city is not None:
            place_shapes = gpd.read_file(texas_places_path)
            places_shapes.crs = {'init':'epsg:4326'}
            city_shape = place_shapes.loc[place_shapes['NAME'] == city, 'geometry'].values[0]

    else:
        result = pd.DataFrame(columns = gdf[factor].unique().tolist())

    if categorical == False and factor is not None:
        print('ignoring nas')
        # To make this play nicely with regulation data, just in case. Can't harm either way.
        gdf = gdf.loc[gdf[factor].notnull()]

    gdf.reset_index(drop = True, inplace = True)
    gdf.index = [str(ind) for ind in gdf.index]
    spatial_index = gdf.sindex


    radius = 0
    while radius <= maximum - step:

        radius += step
        circle = center.buffer(feet_to_mile*radius).difference(center.buffer(feet_to_mile*(radius - step)))

        # Find intersections, and be efficient for larger circles
        if radius <= maximum/2:
            possible_matches_index = list(spatial_index.intersection(circle.bounds))
            possible_matches = gdf.iloc[possible_matches_index]
            # Get precise matches and data
            precise_matches_index = possible_matches.loc[:, geometry_column].intersects(circle)
            precise_matches = possible_matches.loc[precise_matches_index]
            precise_matches.loc[:, geometry_column] = precise_matches.loc[:, geometry_column].intersection(circle) # And adjust geometry
        else:
            grid = fragment(circle, horiz = 10, vert = 10)
            possible_matches_index = set()
            for polygon in grid:
                possible_matches_index = possible_matches_index.union(set(spatial_index.intersection(polygon.bounds)))
            possible_matches = gdf.iloc[list(possible_matches_index)]
            # Get precise matches and data
            precise_matches_index = possible_matches.loc[:, geometry_column].intersects(circle)
            precise_matches = possible_matches.loc[precise_matches_index]
            precise_matches.loc[:, geometry_column] = precise_matches.loc[:, geometry_column].intersection(circle) # And adjust geometry

        precise_matches.loc[:, 'area'] = precise_matches.loc[:, geometry_column].area
        total_area = precise_matches['area'].sum()

        # In this case, the denominator is the total area within the ring
        if factor is None:
            if city is not None:
                area = circle.intersection(city_shape).area
            else:
                area = circle.area
            result.loc[radius] = total_area/area

        # In all other cases, the denominator is the total area inside the 'factor' column
        elif categorical:
            result.loc[radius] = precise_matches[['area', factor]].groupby(factor)['area'].sum().divide(total_area)
        else:
            result.loc[radius] = precise_matches['area'].dot(precise_matches[factor])/(total_area)

    return result

def points_intersect_polygon(points, polygon, spatial_index, points_geometry_column = 'geometry', factor = None, categorical = True, by = 'mean'):
    """
    Given many points and a polygon, finds one of three things. (1) If factor = None, the number of points inside the
    polygon, (2) if factor is not None and categorical = True, the number of points inside the polygon subsetted by a
    categorical factor, (3) if factor is not None and categorical = False, the summarized value (mean/median)
    of a factor associated with each point of each point inside the polygon.
    :param points: A GDF with a geometry column
    :param polygon: The polygon to see whether the points are inside.
    :param spatial_index: The spatial index of the points
    :param factor: The factor to average over or subset by (if categorical).
    :param categorical: If True, then the factor should be treated as a categorical variable.
    :param by: If categorical is False, can either summarize using by = 'mean' or by = 'median'
    :return: float or pandas series
    """
    # Get intersections
    possible_matches_index = list(spatial_index.intersection(polygon.bounds))
    possible_matches = points.iloc[possible_matches_index]
    precise_matches_index = possible_matches[points_geometry_column].intersects(polygon)

    if factor is None:
        return sum(precise_matches_index)
    else:
        precise_matches = points.loc[precise_matches_index[precise_matches_index.values].index.tolist()]
        if categorical == True:
            return precise_matches.groupby(factor)[points_geometry_column].count()
        elif by == 'mean':
            return precise_matches.groupby(factor)[points_geometry_column].mean()
        elif by == 'median':
            return precise_matches.groupby(factor)[points_geometry_column].median()
        else:
            print(
            'In point_choropleth call, "by" argument must either equal "mean" or "median" - you put "{}"'.format(by))
            return None





if __name__ == '__main__':

    #from helpers import process_zoning_shapefile

    #block_data = get_block_geodata(['X01_AGE_AND_SEX'], cities = 'Austin')
    #austin_zones = process_zoning_shapefile(austin_inputs)

    #q = get_all_averages_by_area(block_data, austin_zones, fillna = 0)
    #q['dens'] = q['B01001e1'].divide(q['geometry'].area)
    #q.plot(column = 'B01001e1', legend = True)
    #plt.show()




    import sys
    sys.exit()


    # For testing simplifying function -------------------------------------------------


    # Get initial dataset in several steps

    # Step 1: Get centroids and index
    #data = data[['broad_zone', 'geometry']]
    spatial_index = data.sindex

    # Step 2: Find test polygon and only consider data closest to it
    n = 208
    test_polygon = data.ix[n, 'geometry']
    subset = list(spatial_index.nearest(test_polygon.bounds, num_results=1500))  # .16 milliseconds on average
    data = data.iloc[subset]

    # Step 3: Reset data's index for convenience)
    data.index = np.arange(0, len(data.index), 1)
    data.plot(column = 'broad_zone', legend = True)

    # Run
    data1 = combine_all(data, 'broad_zone', max_comb = 10, alpha = 0.000001, use_v2=True)
    data1.set_geometry('geometry').plot(column = 'broad_zone', legend = True, alpha = 0.5)

    result = process_combined_result(data1, factor = 'broad_zone')
    result.set_geometry('geometry').plot(column = 'broad_zone', legend = True, alpha = 0.5)


    #data1 = combine_all(data1, 'broad_zone', max_comb = 5, alpha = 0.000001, use_v2=True)
    #data1.set_geometry('geometry').plot(column = 'broad_zone', legend = True, alpha = 0.5)

    #data2 = combine_all(data1, 'broad_zone', max_comb = 25, alpha = 400, use_v2=True)
    #data2.set_geometry('geometry').plot(column = 'broad_zone', legend = True, alpha = 0.5)
    #data3 = combine_all(data, 'broad_zone', max_comb = 25, alpha = 0.5)
    #data3.set_geometry('geometry').plot(column = 'broad_zone', legend = True, alpha = 0.5)


    plt.show()





