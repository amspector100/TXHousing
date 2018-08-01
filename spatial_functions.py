# -----------------------------------------------------------------------------------------------------------------
# This module broadly contains a variety of spatial helper functions for the TX Housing project. It has six subsections.
# (1) Some very basic spatial processing functions (get rid of invalid or null geometries, for example)
# (2) Slightly more complex spatial utilities (fragment large polygons, create grids out of points)
# (3) Functions which find distances/areas for entire gdfs without employing spatial joins
# (4) Functions which deal with block geodata and joining it (spatially) with non-block data
# (5) Functions used to create graphs that group polygons/points by their distance from the city center.
# (6) Broader functions used to calculate intersections very quickly. (Every function in this module employees rtrees
#     for intersections, but functions in this section are slightly broader).
# -----------------------------------------------------------------------------------------------------------------

import time
from tqdm import tqdm

import scipy
import shapely
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from inputs import *
import helpers

# Part 1: Processing functions ----------------------------------------------------------------------------------------

# For polygon gdfs: ignore multipolygons and make invalid polygons valid
def process_geometry(gdf, geometry_column = 'geometry', drop_multipolygons = True):
    gdf.loc[:, geometry_column] = gdf[geometry_column].apply(lambda poly: poly if poly.is_valid else poly.buffer(0))
    if drop_multipolygons:
        gdf = gdf.loc[gdf[geometry_column].apply(lambda x: isinstance(x, shapely.geometry.polygon.Polygon))]
    return gdf

# For point gdfs: drop invalid points (e.g., points with Na's instead of coordinates)
def process_points(points, geometry_column = 'geometry'):
    points = points.loc[points[geometry_column].is_valid]
    points.reset_index(drop=True, inplace=True)
    points.index = [str(ind) for ind in points.index]
    return points

# Part 2: Spatial utilities -------------------------------------------------------------------------------------------

# Helper function for spatial tree efficiency - fragment polygon into smaller sizes
def fragment(polygon, horiz = 10, vert = 10):
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

# Get urban cores (i.e. circle of certain radius around a lat long)
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

# Given a bunch of points, make_point_grid splits them into a rectangular grid of m by n and calculates either the
# number of points in each rectangle or the mean or median of a factor associated with the points for each rectangle.

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

 # Part 3: Functions for finding distances/areas for entire gdfs w/out using spatial joins ----------------------------

# Calculates distance between two lat/long points using haversine formula.
# See https://gis.stackexchange.com/questions/279109/calculate-distance-between-a-coordinate-and-a-county-in-geopandas
def haversine(point1, point2, lon1 = None, lat1 = None, lon2 = None, lat2 = None):
    """
    :param point1: Shapely point. Long then lat.
    :param point2: Shapely point. Long then lat.
    :param lon1, lat1, lon2, lat2: Alternatively, supply the long and lattitudes.
    :return: Distance in miles.
    """

    # Retrieve lat and long
    if lon1 is None or lat1 is None:
        lon1, lat1 = list(point1.coords[:][0][0:])
    if lon2 is None or lat2 is None:
        lon2, lat2 = list(point2.coords[:][0][0:])

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(scipy.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = scipy.sin(dlat / 2) ** 2 + scipy.cos(lat1) * scipy.cos(lat2) * scipy.sin(dlon / 2) ** 2
    c = 2 * scipy.arcsin(scipy.sqrt(a))
    r = 3956  # Radius of earth in miles. Use 6371 for km
    return c * r

def calculate_dist_to_center(gdf, zoning_input, drop_centroids = True):
    """
    :param gdf: A GeoDataFrame, in lat/long
    :param zoning_input: Zoning input, used to find the center of the city
    :param drop_centroids: Boolean, default true. If true, drop the centroids inplace after caclulation.
    :return: Pandas Series of floats (distances from center).
    """

    center = shapely.geometry.point.Point(zoning_input.long, zoning_input.lat)

    def dist_to_center(point):
        dist = haversine(point, center)
        return dist

    if 'centroids' not in gdf.columns:
        gdf['centroids'] = gdf['geometry'].centroid

    distances = gdf['centroids'].apply(dist_to_center)

    if drop_centroids:
        gdf.drop('centroids', inplace=True, axis=1)
    return distances

# Get area in miles.
def get_area_in_units(gdf, geometry_column = 'geometry', newproj = 'epsg:2277', scale = 3.58701*10**(-8), name = 'area', final_projection = None):
    """
    Get area of polygons of a geodataframe in units. By default, gets it in miles.
    :param gdf: Geodataframe with polygons in the geometry column.
    :param geometry_column: Geometry column of the geodataframe, defaults to 'geometry'
    :param newproj: The new projection to use to calculate units. Defaults to epsg:2277, which is probably fine for
    Austin/Dallas/Houston and is in feet.
    :param scale: A scale to multiply by. Defaults to 3.58701*10**(-8) which is the number of square miles in a square
    foot.
    :param name: the name of the new column that will be created to store the area information. Defaults to 'area'.
    :param final_projection: The final projection that the returned gdf should be in. Defaults to the gdf's current crs.
    :return: The geodataframe with a column named name (defualts to 'area') which has the area of each polygon in
    the desired units.
    """
    if final_projection is None:
        final_projection = gdf.crs
    gdf = gdf.to_crs({'init':newproj})
    gdf[name] = scale*gdf[geometry_column].area
    gdf = gdf.to_crs(final_projection)
    return gdf

# Part 4: Functions which deal with block data and spatial joins with non-block data ----------------------------------

# See datadic at https://www2.census.gov/geo/tiger/TIGER_DP/2016ACS/Metadata/BG_METADATA_2016.txt
def get_block_geodata(data_layers, cities = None, get_percent_residential = True):
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
        print("""Note: in get_all_averages_by_area, the crs for data_source ({}) and other_geometries ({}) disagree.
              This might cause issues.""".format(data_source.crs, other_geometries.crs))

    # Process data for convenience (just to prevent multipolygons/invalid polygons from messing things up)
    data_source = process_geometry(data_source, drop_multipolygons=drop_multipolygons)
    data_source.reset_index(drop = True)
    data_source.index = [str(ind) for ind in data_source.index]

    # Get the feature in terms of units per area (do a bit of renaming to make it clear these are densities)
    old_columns_dictionary = {str(feature) + '_density':feature for feature in features}
    new_columns_dictionary = {feature:str(feature) + '_density' for feature in features}
    new_columns = [new_columns_dictionary[key] for key in new_columns_dictionary]

    # Defaults to percent residential
    if account_method == 'percent_residential':
        print('Using percent_residential accounting method')
        data_source = data_source.loc[data_source['percent_residential'] >= 0.01] # Ignore blocks which are < .5% resid
        densities = data_source[features].divide(data_source['percent_residential'], axis = 0)
        densities = densities.divide(data_source[data_source_geometry_column].area, axis = 0)
    elif account_method == 'water':
        print('Using water accounting method')
        data_source.loc[:, 'percent_land'] = data_source['ALAND'].divide(data_source["ALAND"] + data_source['AWATER'])
        densities = data_source[features].multiply(data_source['percent_land'], axis = 0)
        densities = densities.divide(data_source[data_source_geometry_column].area, axis = 0)
    else:
        print('Not using any accounting method')
        densities = data_source[features].divide(data_source[data_source_geometry_column].area, axis = 0)

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
        result = get_average_by_area(data_source, spatial_index, polygon, density_features = new_columns,
                                     geometry_column = data_source_geometry_column)
        if fillna is not None:
            result.fillna(fillna, inplace = True)
        return result

    # Get averages by area - this takes a while.
    final_values = other_geometries[other_geometries_column].apply(get_avg)
    final_values = final_values.rename(columns = old_columns_dictionary)
    other_geometries = other_geometries.join(final_values)
    print("Finished calculating average {} by area, took {}".format(features, time.time() - time0))

    return other_geometries

# Part 5: Create dist from city center graphs, for points and polygons ------------------------------------------------

def points_intersect_rings(gdf, zoning_input, factor = None, step = 1, categorical = True, by = 'mean',
                           geometry_column = 'geometry', per_square_mile = True, maximum = None):

    """
    :param gdf: Geodataframe. Assumes this is already transformed to lat long coord system. Should be polygon geometry.
    :param factor: A factor to calculate for. If none, will just calculate the area or number of points radiating out
    from city center.
    :param zoning_input: Zoning input, at the center of which we start our rings
    :param step: Number of miles where the ring radiates outwards.
    :param maximum: Max radius (miles)
    :param categorical: If true, will only calculate percent land (or % of points) in the radius
    :param by: Defaults to "mean". If categorical = False, use "by" to determine how to calculate averages over points.
    :param geometry_column: name of geometry column. Default geometry.
    :param psm: if true, divide by the area of the ring.
    :param maximum: float, defualts to None. If not None, will group everything greater than this maximum into a single
    category.
    :return: Dataframe
    """

    # Get distance from center of city, in miles currently
    center = shapely.geometry.point.Point(zoning_input.long, zoning_input.lat)

    def rounded_dist_to_center(point):
        dist = haversine(point, center)
        rdist = step*((dist // step) + 1) # Do this to get smoother bins (i.e. always round up)
        return rdist

    gdf['dist_to_center'] = gdf[geometry_column].apply(rounded_dist_to_center)

    # Group ouotliers together if told
    if maximum is not None:
        def group_outliers(dist_to_center):
            if dist_to_center > maximum:
                return str(maximum) + '+'
            else:
                return dist_to_center
        gdf['dist_to_center'] = gdf['dist_to_center'].apply(group_outliers)
        if per_square_mile:
            Warning('In points_intersect_rings, dropping all points about the maximum because per_square_mile is True.')
            gdf = gdf.loc[gdf['dist_to_center'].apply(helpers.will_it_float)]


    # Get counts if no factor is provided
    if factor is None:
        result = gdf[['geometry', 'dist_to_center']].groupby('dist_to_center')['geometry'].count()

    # Get counts by factor level if factor is categorical
    elif categorical:
        result = gdf[['dist_to_center', factor]].groupby(['dist_to_center', factor]).size()
        result = result.unstack().fillna(0)

    # Get mean and median else
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
            return math.pi*(radius**2 - (radius - step)**2)

        areas = pd.Series(result.index.map(get_area_ring).tolist(), index = result.index)

        if isinstance(result, pd.DataFrame):
            result = result.divide(areas, axis = 0)
        else:
            result = result.divide(areas)

    # Return result
    return result


# Radius from city center for polygon data
def polygons_intersect_rings(gdf, zoning_input, factor = None, newproj = 'epsg:2277', step = 1,
                             maximum = 20, categorical = True, geometry_column = 'geometry',
                             group_outliers = True, outlier_maximum = 35, city = None):

    """
    :param gdf: Geodataframe. Assumes this is already transformed to lat long coord system. Should be polygon geometry.
    :param factor: A factor to calculate for. If none, will just calculate the area or number of points radiating out
    from city center.
    :param zoning_input: Zoning input, at the center of which we start our rings
    :param newproj: the new projection system necessarily used in this. Defaults to 2277 which is good for Austin and
    fine for Texas. Note units in this are in feet.
    :param step: Number of miles where the ring radiates outwards.
    :param maximum: Max radius (miles)
    :param points: Boolean. If true, assumes geometry of gdf is a point geometry.
    :param geometry_column: name of geometry column. Default geometry.
    :param categorical: If true, will only calculate percent land (or % of points) in the radius. Else will calculate
    mean by area.
    :param city: If not "none", if factor is "none", will read the shaepfile of the boundaries of the city to ensure
    more accurate calculations. (Otherwise, for a radius of size 12, the area of the circle might be greater than the
    area of the city).
    :param group_outliers: Boolean, defaults to true. If true, group everything with a distance greater than the maximum
     into one group (of maximum size).
    :param outlier_maximum: Float, defaults to 35. For computational efficiency, this function will not consider outliers
     higher than this distance from the cneter of the city.
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

    # If this optional arg is true, then lump everything else into one category
    if group_outliers:
        circle = center.buffer(feet_to_mile*outlier_maximum).difference(center.buffer(feet_to_mile*radius))
        grid = fragment(circle, horiz=15, vert=15)
        possible_matches_index = set()
        for polygon in grid:
            possible_matches_index = possible_matches_index.union(set(spatial_index.intersection(polygon.bounds)))
        possible_matches = gdf.iloc[list(possible_matches_index)]
        # Get precise matches and data
        precise_matches_index = possible_matches.loc[:, geometry_column].intersects(circle)
        precise_matches = possible_matches.loc[precise_matches_index]
        precise_matches.loc[:, geometry_column] = precise_matches.loc[:, geometry_column].intersection(circle)  # And adjust geometry

        # Compute areas
        precise_matches.loc[:, 'area'] = precise_matches.loc[:, geometry_column].area
        total_area = precise_matches['area'].sum()

        # In this case, the denominator is the total area within the ring
        label = str(radius) + '+'
        if factor is None:
            if city is not None:
                area = circle.intersection(city_shape).area
            else:
                area = circle.area
            result.loc[label] = total_area/area

        # In all other cases, the denominator is the total area inside the 'factor' column
        elif categorical:
            result.loc[label] = precise_matches[['area', factor]].groupby(factor)['area'].sum().divide(total_area)
        else:
            result.loc[label] = precise_matches['area'].dot(precise_matches[factor])/(total_area)

    return result

# Part 6: Fast intersection functions and some slight variations on them ----------------------------------------------

def points_intersect_polygon(points, polygon, spatial_index, points_geometry_column = 'geometry',
                             factor = None, categorical = True, by = 'mean'):
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

def fast_polygon_intersection(gdf, large_polygon_list, geometry_column = 'geometry', names = None, **kwargs):
    """
    # A very simple function here - quickly determines which of a set of small polygons, which should be part of a gdf,
    are inside a list of large polygons (this should be a list or a GeoSeries). Uses centroids to speed up computation.
    :param large_polygon_list: List of large polygons (i.e. municipal boundaries)
    :param names: A list of the names of the large polygons, in the same order as the large_polygon_list.
    :param gdf: A geodataframe with polygon geometries.
    :param geometry_column: The geometry column of the geodataframe.
    :param **kwargs: kwargs to pass to the "fragment" call (fragmenting large polygons into grids speeds up computation
    due to the nature of rtrees).
    :return: If names = None, a list of the indexes of the gdf of the small polygons which lie inside at least one of
    the large polygons. Otherwise, will return a dictionary which maps the indexes of the gdp of small polygons to
    the names of the large polygons. Order matters here: if a point lies in multiply large polygons (which is not what
    this function is intended for), the dictionary will map the index of the point to the name of the last
    """

    # Convert geoseries input to list
    if isinstance(large_polygon_list, gpd.GeoSeries):
        large_polygon_list = large_polygon_list.values.tolist()

    # Get centroids
    gdf['centroids'] = gdf[geometry_column].centroid
    gdf = gdf.set_geometry('centroids')
    spatial_index = gdf.sindex

    # Initialize result
    result = {}
    if names is None:
        names_list = [0]*len(large_polygon_list)
    else:
        names_list = names

    # Loop through and find intersections
    for name, polygon in zip(names_list, large_polygon_list):
        grid = fragment(polygon, **kwargs)
        for grid_piece in grid:
            possible_intersections_index = list(spatial_index.intersection(grid_piece.bounds))
            possible_intersections = gdf.iloc[possible_intersections_index]
            precise_intersections_bools = possible_intersections['centroids'].intersects(grid_piece)
            precise_intersections = possible_intersections[precise_intersections_bools].index.tolist()
            for i in precise_intersections:
                result[i] = name

    # If the names are meaningless, just return the list of the result
    if names is None:
        result = [key for key in result]

    return result

# Intersect zoning with zip codes with great precision
def zip_intersect(gdf, zips):
    """
    :param gdf: A geodataframe of some sort (it should probably be in the US, otherwise this is pointless).
    :param zips: The zip codes in the area of interest.
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

# Returns the number of points inside the polygons divided by the area of the polygons in square miles
def points_over_area(points, polygons, points_geometry_column='geometry', poly_geometry_column='geometry',
                     normalize_by_area=True):
    """
    :param points: GeoDataframe of points, should be in lat long
    :param polygons: Geodataframe of polygons
    :param: normalize_by_area: If true, will divide the number of points by the area of the location (in square miles).
    :return: GeoDataFrame of polygons
    """

    # Process points and polys to make them valid
    polygons = process_geometry(polygons, drop_multipolygons=False)
    polygons.reset_index(drop=True, inplace=True)
    polygons.index = [str(ind) for ind in polygons.index]

    points = process_points(points, geometry_column=points_geometry_column)

    counter = 0

    # Intersections (efficiently)
    spatial_index = points.sindex
    for poly in polygons[poly_geometry_column]:
        possible_intersections_index = list(spatial_index.intersection(poly.bounds))
        possible_intersections = points.iloc[possible_intersections_index]
        precise_matches = possible_intersections[points_geometry_column].intersects(poly)
        number_matches = sum(precise_matches.tolist())
        counter += number_matches

    # Return # of points divided by area
    area = get_area_in_units(polygons)['area'].sum()

    if normalize_by_area:
        return counter / area  # Would be nice to get this area in miles. Oh well.
    else:
        return counter
