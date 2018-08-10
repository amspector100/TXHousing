""" This module is meant to simplify very complex parcel and zoning shapefiles to allow graphing and improve rendering.
This is currently NOT used in any other module because it reduces the accuracy of zoning maps.
Currently, they are NOT used in the mastermaps because they take a while to run and substantially reduce the
accuracy of the initial data. It's also untested post-reorganization."""

import math
import time
import numpy as np
import pandas as pd
import geopandas as gpd
import tqdm

import shapely.geometry
from shapely.ops import cascaded_union, polygonize, unary_union
from scipy.spatial import Delaunay

# Global granularity variable
max_comb = 30
alpha = 0.01 # I got this by experimenting

# Adapted from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
def alpha_shape(points, alpha = 0.001):
    """
    Compute the alpha shape/concave hull of points.
    :param points: Iterable container of shapely points.
    :param alpha: alpha value to influence the the border. Smaller numbers
        don't fall inward as much as larger numbers.
    """

    # This algorithm is pointless unless there are more than 3 points
    if len(points) < 4:
        return shapely.geometry.multipoint.MultiPoint(list(points)).convex_hull

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
    m = shapely.geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    result = cascaded_union(triangles)
    if result.is_valid == False:
        result = result.buffer(0)
    return result # could return edge_points

# Combine data
def combine(data, feature, spatial_index, k, max_comb = max_comb, alpha = alpha):

    # Get polygon and zone
    test_polygon = data.loc[k, 'geometry']
    test_zone = data.loc[k, feature]

    # Create index and find polygons
    new_index = spatial_index.nearest(test_polygon.bounds, num_results = max_comb) # .16 milliseconds on average
    polys = data.loc[list(new_index)]
    good_polys = [shapely.geometry.point.Point(coord) for coord in list(test_polygon.exterior.coords)]
    good_inds = [k]

    # Loop through to see which ones are in the right zone
    for ind, row in polys.iterrows():
        if row[feature] == test_zone and row['flag'] == False:
            good_polys.extend([shapely.geometry.point.Point(coord) for coord in list(row['geometry'].exterior.coords)])
            good_inds.append(ind)
        else:
            break

    # Get convex hull and return
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
        line = shapely.geometry.linestring.LineString([test_polygon.centroid.coords[:][0], possible_poly.centroid.coords[:][0]])
        if any([line.intersects(bad_poly) for bad_poly in bad_polys['geometry']]):
            continue
        else:
            good_poly_coords.extend([shapely.geometry.point.Point(coord) for coord in list(row['geometry'].exterior.coords)])
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
            cached_geometry = result.loc[zone, 'geometry'].copy()
            result.at[zone, 'geometry'] = unary_union(cached_geometry.difference(running_union))
            running_union = cached_geometry.union(running_union)

    result[factor] = result.index
    return result