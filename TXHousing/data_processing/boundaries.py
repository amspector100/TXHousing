"""Paths and methods for boundary files (uas, csa, counties, blocks, etc)"""
from ..utilities import measurements, spatial_joins
import geopandas as gpd
import warnings
import shapely.geometry

csa_path = "data/cb_2017_us_csa_500k/cb_2017_us_csa_500k.shp"
cbsa_path =  "data/cb_2017_us_cbsa_500k/cb_2017_us_cbsa_500k.shp"
ua_path = "data/cb_2017_us_ua10_500k/cb_2017_us_ua10_500k.shp"

zip_boundaries_path = "data/cb_2017_us_zcta510_500k/cb_2017_us_zcta510_500k.shp"
county_boundaries_path = "data/cb_2017_us_county_500k/cb_2017_us_county_500k.shp"
texas_blocks_path = "data/ACS_2016_5YR_BG_48_TEXAS.gdb"
texas_places_path = measurements.texas_places_path


class Boundaries():
    """A parent class for boundary files (i.e. counties, blocks, places).

    :param index_col: The column in the shapefile to use as the index."""


    def fast_intersection(self, gdf, geometry_column = 'geometry', **kwargs):
        """

        Given a gdf of polygons or points, will calculate which boundary each polygon/point lies inside. Assumes that
        if the gdf is full of polygon data, each polygon will only intersect a single zip code.

        :param gdf: A geodataframe, in the same crs as the zip data.
        :param geometry_column: The geometry column of the gdf.
        :param **kwargs: kwargs to pass to the underlying fast_polygon_intersection or points_intersect_multiple_polygons
            functions.
        :return: A pandas series mapping the index of the gdf to the names of the zip codes. If an index does not
        appear in the returned series, that is because the point/polygon for that index did not lie inside any of
        the zip codes.

        """

        if gdf.crs is None:
            warnings.warn('In Boundaries.fast_intersection, assuming external gdf data is in lat/long because no crs is given')
        elif gdf.crs != {'init':'epsg:4326'}:
            warnings.warn('In Boundaries.fast_intersection, forced to transform external gdf data to lat/long')

        sample_geometry = gdf.loc[gdf.index[0], geometry_column]
        if isinstance(sample_geometry, shapely.geometry.polygon.Polygon):
            result = spatial_joins.fast_polygon_intersection(small_polygon_gdf = gdf, large_polygon_gdf = self.data,
                                                             small_geometry_column = geometry_column, **kwargs)
            return result
        elif isinstance(sample_geometry, shapely.geometry.point.Point):
            result = spatial_joins.points_intersect_multiple_polygons(points_gdf = gdf, polygon_gdf = self.data,
                                                                      points_geometry_column = geometry_column, **kwargs)
            return result
        else:
            raise TypeError('gdf geometry must be a Shapely point or polygon, not {}'.format(type(sample_geometry)))

    def spatial_summary(self, gdf, geometry_column = 'geometry', **kwargs):

        """

        Given a gdf of polygons or points, will calculate summaries of

        :param gdf: A geodataframe, in the same crs as the zip data.
        :param geometry_column: The geometry column of the gdf.
        :param **kwargs: kwargs to pass to the underlying fast_polygon_intersection or points_intersect_multiple_polygons
            functions.
        :return: A pandas series mapping the index of the gdf to the names of the zip codes. If an index does not
        appear in the returned series, that is because the point/polygon for that index did not lie inside any of
        the zip codes.

        """

# For later use
def __init__(self, ziplist=None):
    self.ziplist = list(map(str, ziplist))

    # Read in data and subset
    self.data = gpd.read_file(zip_boundaries_path)
    self.data.index = self.data['ZCTA5CE10']
    self.data.index = self.data.index.map(str)
    if self.ziplist is not None:
        self.data = self.data.loc[ziplist]






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