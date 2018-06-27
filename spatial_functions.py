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
from shapely import geometry
from shapely.ops import cascaded_union, polygonize, unary_union
from scipy.spatial import Delaunay
import math
from tqdm import tqdm
import helpers
import copy

# Everything in this file should be sufficiently general as to work without modifications for any city in the United States

# Basically meant to simplify very complex zoning data for categorical features ---------------------------------------

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
def process_geometry(gdf):
    gdf.at[:, 'geometry'] = gdf['geometry'].apply(make_valid_buffer)
    gdf = gdf.loc[gdf['geometry'].apply(is_polygon)]
    return gdf

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

# Intersect zoning with zip codes -------------------------------------------------------------------------------------
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

if __name__ == '__main__':

    # For testing in general --------------------------------------------------------------

    from inputs import *
    from helpers import process_zoning_shapefile

    # The data is indexed by a range index (start = 0, stop = 21623)
    data = process_zoning_shapefile(austin_inputs, broaden = True)
    data['broad_zone'] = data['broad_zone'].apply(str)

    # Process data by only considering the polygons - might fix this later, it only excludes 40/21.6K zones though.
    data = process_geometry(data)
    data = data.loc[data['broad_zone'].notnull()]


    # For testing intersection function ---------------------------------------------------
    #from inputs import austin_zips
    #zip_intersect(data, austin_zips)

    #import sys
    #sys.exit()

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





