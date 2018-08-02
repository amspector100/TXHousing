"""Fast intersection functions"""

import pandas as pd
import geopandas as gpd
import warnings
from . import simple

def points_intersect_single_polygon(points, polygon, spatial_index, points_geometry_column = 'geometry',
                             factor = None, categorical = True, by = 'mean', **kwargs):
    """

    Given many points and a polygon, finds one of three things. (1) If factor = None, the number of points inside the
    polygon, (2) if factor is not None and categorical = True, the number of points inside the polygon subsetted by a
    categorical factor, (3) if factor is not None and categorical = False, the summarized value (mean/median)
    of a factor associated with each point of each point inside the polygon.

    :param points: A GDF with a points geometry column
    :param polygon: The polygon to see whether the points are inside.
    :param spatial_index: The spatial index of the points
    :param factor: The factor to average over (if continuous) or subset by (if categorical).
    :param categorical: If True, then the factor should be treated as a categorical variable.
    :param by: If categorical is False, can either summarize using by = 'mean' or by = 'median'
    :param **kwargs: Kwargs to pass to the "fragment" function in the TXHousing.utilities.simple module. Fragmenting polygons
        speeds up the computation for all but very small polygons. If you do not want to fragment the polygons (the
        only reason to do this is speed, it will not affect the results), pass in horiz = 1 and vert = 1 as kwargs.
    :return: float or pandas series

    Note: it is often useful to apply this function to an entire gdf of polygons.

    """

    # Get intersections
    all_precise_matches_indexes = set()
    grid = simple.fragment(polygon, **kwargs)
    for grid_piece in grid:
        possible_matches_index = list(spatial_index.intersection(grid_piece.bounds))
        possible_matches = points.iloc[possible_matches_index]

        # Prevent weird indexing errors by ignoring empty results
        if possible_matches.shape[0] == 0:
            continue

        precise_matches_index = set(possible_matches.index[possible_matches.intersects(grid_piece)].tolist())
        all_precise_matches_indexes = all_precise_matches_indexes.union(precise_matches_index)

    if factor is None:
        return len(precise_matches_index)
    else:
        precise_matches = points.loc[all_precise_matches_indexes]
        if categorical == True:
            return precise_matches.groupby(factor)[points_geometry_column].count()
        elif by == 'mean':
            return precise_matches[factor].mean()
        elif by == 'median':
            return precise_matches[factor].median()
        else:
            raise ValueError(
            'In points_intersect_polygon call, "by" must either equal "mean" or "median," not "{}"'.format(by))


def polygons_intersect_single_polygon(small_polygons, polygon, spatial_index, geometry_column = 'geometry',
                                      factor = None, categorical = True, account_for_area = True, by = 'mean', **kwargs):
    """

    Given many polygons (i.e. parcels) and a larger polygon (i.e. county boundary), finds one of three things.
    (1) If factor = None, the percent area of the large polygon that is covered by the small polygons
    (2) If factor is not None and categorical = True, the percent area of the large polygon that is covered by the small
    polygons conditional on the factor
    (3) if factor is not None and categorical = False, the summarized value (mean/median) of a factor associated with
    each polygon inside the polygon.

    :param small_polygons: A GDF with a polygon geometry column
    :param polygon: The polygon to see whether the small_polygons are inside.
    :param spatial_index: The spatial index of the small_polygons
    :param factor: The factor to average over (if continuous) or subset by (if categorical).
    :param categorical: If True, then the factor should be treated as a categorical variable.
    :param by: If categorical is False, can either summarize using by = 'mean' or by = 'median'
    :param account_for_area: If categorical is False, by = 'mean', and account_for_area = True, then instead of returning
        the mean of the factor, this will return the divided dot product of the mean and the area of each small_polygon
        that intersects the large_polygon.
    :param **kwargs: Kwargs to pass to the "fragment" function in the TXHousing.utilities.simple module. Fragmenting polygons
        speeds up the computation for all but very small polygons. If you do not want to fragment the polygons (the
        only reason to do this is speed, it will not affect the results), pass in horiz = 1 and vert = 1 as kwargs.
    :return: float (if categorical = False or factor is None) or pandas series

    Note: it is often useful to apply this function to an entire gdf of polygons.

    """

    # Get the indexes of the small polygons which intersect the big polygon
    precise_matches_indexes = set()
    grid = simple.fragment(polygon, **kwargs)
    for grid_piece in grid:
        possible_matches_index = list(spatial_index.intersection(polygon.bounds))
        possible_matches = small_polygons.iloc[possible_matches_index]

        # Might speed things up a tiny bit; also might prevent weird indexing errors
        if possible_matches.shape[0] == 0:
            continue

        precise_matches_to_add = set(possible_matches.index[possible_matches.intersects(grid_piece)].tolist())
        precise_matches_indexes = precise_matches_indexes.union(precise_matches_to_add)

    # Calculate actual intersections
    precise_matches = small_polygons.loc[precise_matches_indexes]
    precise_matches.loc[:, geometry_column] = precise_matches.loc[:, geometry_column].intersection(polygon)
    precise_matches['area'] = precise_matches[geometry_column].area

    # Calculate answers
    if factor is None:
        return precise_matches['area'].sum()
    else:
        if categorical == True:
            return precise_matches.groupby(factor)['area'].sum().divide(polygon.area)
        elif by == 'mean':
            if account_for_area:
                return precise_matches[factor].multiply(precise_matches['area']).sum()/(polygon.area)
            else:
                return precise_matches[factor].mean()
        elif by == 'median':
            return precise_matches[factor].median()
        else:
            raise ValueError(
            'In points_intersect_polygon call, "by" must either equal "mean" or "median," not "{}"'.format(by))


def points_intersect_multiple_polygons(points_gdf, polygons_gdf, points_spatial_index = None,
                              points_geometry_column = 'geometry', polygons_geometry_column = 'geometry',
                              polygons_names_column = None, **kwargs):
    """

    Given a gdf of points and a gdf of polygons, calculates the polygon in which each point lies. This function assumes
    that each point will lie in at most one of the polygons. If that assumption is not true, use instead the
    points_intersect_single_polygon function and apply it to the geometry column of a polygon gdf.

    :param points_gdf: A geodataframe of points data.
    :param polygons_gdf:
    :param points_spatial_index: Optional; the spatial_index of the points geodataframe. If not supplied, the function
        will automatically generate the spatial index.
    :param points_geometry_column: Geometry column for the points data.
    :param polygons_geometry_column: Geometry column for the polygon data.
    :param polygons_names_column: Column for the names of each polygon; if none will use the index of the polygons_gdf.
    :param kwargs: Kwargs to pass to the "fragment" function in the TXHousing.utilities.simple module. Fragmenting polygons
        speeds up the computation for all but very small polygons. If you do not want to fragment the polygons (the
        only reason to do this is speed, it will not affect the results), pass in horiz = 1 and vert = 1 as kwargs.
    :return: A pandas series mapping the index of the points_gdf to the names of the polygons. If an index does not
        appear in the returned series, that is because the point corresponding to that index did not lie inside any of
        the polygons.
    """

    result = {}
    if points_spatial_index is None:
        points_spatial_index = points_gdf.sindex

    # Create lists of polygons/names from large_polygon_gdf
    polygon_list = polygons_gdf[polygons_geometry_column].values.tolist()
    if polygons_names_column is not None:
        names_list = polygons_gdf[polygons_names_column].values.tolist()
    else:
        names_list = polygons_gdf.index.tolist()

    # Loop through and find intersections
    warning_count = 0
    for name, polygon in zip(names_list, polygon_list):

        # Fragment the polygons.
        grid = simple.fragment(polygon, **kwargs)
        for grid_piece in grid:
            possible_intersections_index = list(points_spatial_index.intersection(grid_piece.bounds))
            possible_intersections = points_gdf.iloc[possible_intersections_index]
            precise_intersections_bools = possible_intersections[points_geometry_column].intersects(grid_piece)
            precise_intersections = possible_intersections[precise_intersections_bools].index.tolist()
            for i in precise_intersections:
                if i in result and result[i] != name:
                    warning_count += 1
                result[i] = name

    # Warn the user if a point lies in multiple polygons
    if warning_count != 0:
        warnings.warn("""In points_intersect_polygons, up to {} points are in multiple polygons, but points_intersect_polygons
         only returns one polygon per point (if a point is in two polygons, it will only show up as being in one).""".format(warning_count))

    result = pd.Series(result)
    return result


def fast_polygon_intersection(small_polygon_gdf, large_polygon_gdf, small_geometry_column = 'geometry', large_geometry_column = 'geometry',
                              large_name_column = None, **kwargs):
    """

    Given a gdf of small polygons (i.e. parcels) and a gdf of large polygons (i.e. municipal boundaries), calculates the
    large polygon in which each small polygon lies. This function is based on the points_intersect_multiple_polygons
    function and therefore assumes that each small polygon will lie in at most one of the large polygons.

    :param small_polygon_gdf: A gdf of small polygons (i.e. parcels)
    :param large_polygon_gdf: A gdf of large polygons (i.e. municipal boundaries)
    :param small_geometry_column: The geometry column of the small_polygon_gdf.
    :param large_geometry_column: The geometry column of the large_polygon_gdf.
    :param large_name_column: Column for the names of each large polygon; if none will use the index of the large_polygon_gdf.
    :param kwargs: Kwargs to pass to the "fragment" function in the TXHousing.utilities.simple module. Fragmenting polygons
        speeds up the computation for all but very small polygons. If you do not want to fragment the polygons (the
        only reason to do this is speed, it will not affect the results), pass in horiz = 1 and vert = 1 as kwargs.
    :return: A pandas series mapping the index of the small polygons to the names of the large polygons. If an index
        does not appear in the returned series, that is because the small polygon corresponding to that index did not
        lie inside any of the large polygons.
    """

    # Get centroids
    small_polygon_gdf['centroids'] = small_polygon_gdf[small_geometry_column].centroid
    small_polygon_gdf = small_polygon_gdf.set_geometry('centroids')

    result = points_intersect_multiple_polygons(points_gdf = small_polygon_gdf, polygons_gdf = large_polygon_gdf,
                                                points_geometry_column = 'centroids',
                                                polygons_geometry_column = large_geometry_column,
                                                polygons_names_column = large_name_column, **kwargs)

    small_polygon_gdf.set_geometry(small_geometry_column) # Undo global effects on small_polygon_gdf

    return result