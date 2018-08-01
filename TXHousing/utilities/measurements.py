"""Functions which focus on efficienctly measuring distances and area."""

import scipy
import shapely
import pandas as pd
import geopandas as gpd
import simple
import warnings

texas_places_path = "data/cb_2017_48_place_500k/cb_2017_48_place_500k.shp"

# Haversine, area, dist to center - None of these use spatial joins -------------------------------------------------

def haversine(point1, point2, lon1 = None, lat1 = None, lon2 = None, lat2 = None):
    """
    Haversine function calculates distance (in miles) between two points in lat/long coordinates. See
        https://gis.stackexchange.com/questions/279109/calculate-distance-between-a-coordinate-and-a-county-in-geopandas

    :param point1: Shapely point. Long then lat.
    :param point2: Shapely point. Long then lat.
    :param lon1, lat1, lon2, lat2: Alternatively, supply the longitudes and lattitudes manually.
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

def calculate_dist_to_center(gdf, lat, long, drop_centroids = True):
    """

    Calculates distance to the center of the city using haversine on the centroids of objects

    :param gdf: A GeoDataFrame, with lat/long crs. Can either have point or polygon geometry.
    :param lat: Latitude of the center of the city
    :param long: Longitude of the center of the city.
    :param drop_centroids: Boolean, default true. If true, drop the centroids inplace after calculation.
    :return: Pandas Series of floats (distances from center).
    """

    center = shapely.geometry.point.Point(long, lat)

    def dist_to_center(point):
        dist = haversine(point, center)
        return dist

    if 'centroids' not in gdf.columns:
        gdf['centroids'] = gdf['geometry'].centroid

    distances = gdf['centroids'].apply(dist_to_center)

    if drop_centroids:
        gdf.drop('centroids', inplace=True, axis=1)
    return distances

def get_area_in_units(gdf, geometry_column = 'geometry', newproj = 'epsg:2277', scale = 3.58701*10**(-8), name = 'area', final_projection = None):
    """

    Get the area of each polygon of a geodataframe in units of your choice (defaults to square miles). This function
        relies on crs transformations, so for large/complex gdfs, this function is very computationally expensive.

    :param gdf: Geodataframe with polygons in the geometry column.
    :param geometry_column: Geometry column of the geodataframe, defaults to 'geometry'
    :param newproj: The new projection to use to calculate units. Defaults to epsg:2277, which is probably fine for
    Austin/Dallas/Houston and is in feet.
    :param scale: A scale to multiply by. Defaults to 3.58701*10**(-8) which is the number of square miles in a square
    foot.
    :param name: The name of the new column that will be created to store the area information. Defaults to 'area'.
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

# These functions group gdfs into successive rings in distance from city center. They employ spatial joins and are
# somewhat expensive.


# Part 5: Create dist from city center graphs, for points and polygons ------------------------------------------------

def points_intersect_rings(gdf, lat, long, factor = None, step = 1, categorical = True, by = 'mean',
                           geometry_column = 'geometry', per_square_mile = True, maximum = None):

    """
    :param gdf: Geodataframe. Assumes this is already transformed to lat long coord system. Should be polygon geometry.
    :param factor: A factor to calculate for. If none, will just calculate the area or number of points radiating out
    from city center.
    :param lat: The latitude of the center of the rings.
    :param long: The longitude of the center of the rings.
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
    center = shapely.geometry.point.Point(long, lat)

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
            warnings.warn('In points_intersect_rings, dropping all points about the maximum because per_square_mile is True.')
            gdf = gdf.loc[gdf['dist_to_center'].apply(simple.will_it_float)]


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
        raise TypeError("""In points_intersect_rings call, when categorical = False and factors is not None, "by" argument must
        either equal "mean" or "median" - you put "{}" """.format(by))

    # If told, get area per square mile
    if per_square_mile:

        def get_area_ring(radius):
            return scipy.pi*(radius**2 - (radius - step)**2)

        areas = pd.Series(result.index.map(get_area_ring).tolist(), index = result.index)

        if isinstance(result, pd.DataFrame):
            result = result.divide(areas, axis = 0)
        else:
            result = result.divide(areas)

    # Return result
    return result


# Radius from city center for polygon data
def polygons_intersect_rings(gdf, lat, long, factor = None, newproj = 'epsg:2277', step = 1,
                             maximum = 20, categorical = True, geometry_column = 'geometry',
                             group_outliers = True, outlier_maximum = 35, city = None):

    """
    :param gdf: Geodataframe. Assumes this is already transformed to lat long coord system. Should be polygon geometry.
    :param factor: A factor to calculate for. If none, will just calculate the area or number of points radiating out
    from city center.
    :param lat: The latitude of the center of the rings.
    :param long: The longitude of the center of the rings.
    :param newproj: the new projection system necessarily used in this. Defaults to 2277 which is good for Austin and
    fine for Texas. Note units in this are in feet.
    :param step: Number of miles where the ring radiates outwards.
    :param maximum: Max radius (miles)
    :param points: Boolean. If true, assumes geometry of gdf is a point geometry.
    :param geometry_column: name of geometry column. Default geometry.
    :param categorical: If true, will only calculate percent land (or % of points) in the radius. Else will calculate
    mean by area.
    :param city: If not "none", if factor is "none", will read the shapefile of the boundaries of the city to ensure
    more accurate calculations. (Otherwise, for a radius of size 12, the area of the circle might be greater than the
    area of the city).
    :param group_outliers: Boolean, defaults to true. If true, group everything with a distance greater than the maximum
     into one group (of maximum size).
    :param outlier_maximum: Float, defaults to 35. For computational efficiency, this function will not consider outliers
     higher than this distance from the cneter of the city.
    :return: Dataframe or Series
    """

    feet_to_mile = 5280
    gdf = simple.process_geometry(gdf)
    gdf.reset_index(drop = True)
    gdf.index = [str(ind) for ind in gdf.index]


    # Necessarily will have to transform into new epsg code in order to efficiently calculate this.
    # Here, gdf should be in (lat, long) format
    if gdf.crs is None:
        warnings.warn('No CRS set, assuming lat and long data before transformation.')
        gdf.crs = {'init':'epsg:4326'}

    gdf = gdf.to_crs({'init':newproj})


    # Get center and initialize output
    center_latlong = shapely.geometry.point.Point(long, lat)
    center_series = gpd.GeoSeries(center_latlong)
    center_series.crs = {'init':'epsg:4326'}
    center_gdf = gpd.GeoDataFrame(geometry = center_series)
    center_gdf = center_gdf.to_crs({'init':newproj})
    center = center_gdf.loc[0, 'geometry']


    # Initialize result. If categorical, need a dataframe (one column for each unique value). Else, use a pd.Series.
    if factor is None or categorical == False:
        result = pd.Series()
        if city is not None and factor is None:
            place_shapes = gpd.read_file(texas_places_path)
            place_shapes.crs = {'init':'epsg:4326'}
            city_shape = place_shapes.loc[place_shapes['NAME'] == city, 'geometry'].values[0]

    else:
        result = pd.DataFrame(columns = gdf[factor].unique().tolist())

    if categorical == False and factor is not None:
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
            grid = simple.fragment(circle, horiz = 10, vert = 10)
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
        elif total_area == 0:
            result.loc[radius] = float('NaN')

        # In all other cases, the denominator is the total area inside the 'factor' column
        elif categorical:
            result.loc[radius] = precise_matches[['area', factor]].groupby(factor)['area'].sum().divide(total_area)
        else:
            result.loc[radius] = precise_matches['area'].dot(precise_matches[factor]) / (total_area)

    # If this optional arg is true, then lump everything else into one category
    if group_outliers:
        circle = center.buffer(feet_to_mile*outlier_maximum).difference(center.buffer(feet_to_mile*radius))
        grid = simple.fragment(circle, horiz=15, vert=15)
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
        elif total_area == 0:
            result.loc[radius] = float('NaN')

        # In all other cases, the denominator is the total area inside the 'factor' column
        elif categorical:
            result.loc[radius] = precise_matches[['area', factor]].groupby(factor)['area'].sum().divide(total_area)
        else:
            result.loc[radius] = precise_matches['area'].dot(precise_matches[factor]) / (total_area)

    # If this optional arg is true, then lump everything else into one category
    return result
