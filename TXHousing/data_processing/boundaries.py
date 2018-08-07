"""Paths and methods for boundary files (uas, csa, counties, blocks, etc)"""
from ..utilities import simple, measurements, spatial_joins
import geopandas as gpd
import pandas as pd
import warnings
import shapely.geometry
import os

# One exception is that this is needed for utilities, so the places path is listed there.
texas_places_path = measurements.texas_places_path

# All other boundaries paths are listed here.
csa_path = "data/cb_2017_us_csa_500k/cb_2017_us_csa_500k.shp"
cbsa_path =  "data/cb_2017_us_cbsa_500k/cb_2017_us_cbsa_500k.shp"
ua_path = "data/cb_2017_us_ua10_500k/cb_2017_us_ua10_500k.shp"

zip_boundaries_path = "data/cb_2017_us_zcta510_500k/cb_2017_us_zcta510_500k.shp"
county_boundaries_path = "data/cb_2017_us_county_500k/cb_2017_us_county_500k.shp"
texas_blocks_path = "data/ACS_2016_5YR_BG_48_TEXAS.gdb"

# This is the percent of land in each block group that is zoned as residential land. The function below calculates it for
# Austin and Dallas to improve accuracy of some population-based calculations.
bg_percent_residential_path = 'shared_data/bg_percent_residential.csv'


class Boundaries():
    """A parent class for boundary files (i.e. counties, places, uas). Note that block data has a different __init__
    method because block data comes as a geodatabase not a shapefile.

    :param path: The path of the boundaries shapefile.
    :param crs: Default None. If you cannot read the CRS data from the
    :param index_col: The column in the shapefile to use as the index. Defaults to None.
    :param index_type: Change the index to this type, i.e. 'str'
    :param subset_by: Subset by this column. Can be 'index' or another column name.
    :param subset_to: Subset to one of these values.
    :param bounding_counties: Default None. A list of counties. If not None, __init__ will subset the data
        to only include boundaries intersecting these counties.
    :param bounding_polygon: Default None. A shapely polygon. If not None, __init__ will subset the data to only include
        boundaries intersecting the polygon.
    :param to_latlong: Default True. Initially tansform data to lat long."""


    def spatial_subset(self, polygon):
        """ Given a polygon, will return a subset of self.data which only includes boundaries which intersect the polygon."""
        possible_intersections = self.data.iloc[list(self.sindex.intersection(polygon.bounds))]
        result = possible_intersections.loc[possible_intersections['geometry'].intersects(polygon)]
        return result

    def __init__(self, path, crs = None, index_col = None, index_type = None, subset_by = None, subset_to = None,
                 bounding_counties = None, bounding_polygon = None, to_latlong = True):

        # Read data, get sindex, crs
        self.path = path
        self.data = gpd.read_file(path)
        if self.data.crs is None:
            if crs is not None:
                self.data.crs = None
            else:
                warnings.warn('No CRS provided for boundaries data at {}'.format(self.path))
        if to_latlong:
            self.data = self.data.to_crs({'init':'epsg:4326'})

        self.sindex = self.data.sindex

        # Set index and subset
        if index_col is not None:
            self.data.index = self.data[index_col]
        if index_type is not None:
            self.data.index = self.data.index.astype(index_type)
        if subset_by == 'index':
            self.data = self.data.loc[self.data.index.isin(subset_to)]
        elif subset_by is not None:
            self.data = self.data.loc[self.data[subset_by].isin(subset_to)]

        # Spatial subsets
        if bounding_polygon is not None:
            self.data = self.spatial_subset(bounding_polygon)
        if bounding_counties is not None:
            counties = gpd.read_file(county_boundaries_path)
            county_shape = counties.loc[(counties['NAME'].isin(bounding_counties)) & (counties['STATEFP'] == '48')].simplify(tolerance = 0.005).unary_union
            self.data = self.spatial_subset(county_shape)

    def type_of_external_gdf(self, gdf, geometry_column = 'geometry'):
        """Figures out the type of geometry (point or polygon) of a gdf and raises an error if it's mixed"""

        geom_types = set(gdf[geometry_column].geom_type.unique().tolist())
        if len(geom_types) != 1:
            if geom_types != set(['Polygon', 'MultiPolygon']) and geom_types != set(['Point', 'MultiPoint']):
                raise ValueError('Cannot handle mixed geometry types of {}'.format(geom_types))
        if 'Polygon' in geom_types:
            return 'Polygon'
        if 'Point' in geom_types:
            return 'Point'

    def process_external_gdf(self, gdf, geometry_column = 'geometry', **kwargs):
        """
        Processes external gdf's geometry/crs. Passes kwargs to process_geometry or process_points call.
        """

        # Get rid of invalid geometries
        geometry_type = self.type_of_external_gdf(gdf)
        if geometry_type == 'Polygon':
            gdf = simple.process_geometry(gdf, geometry_column = geometry_column, **kwargs)
        elif geometry_type == 'Point':
            gdf = simple.process_points(gdf, geometry_column = geometry_column, **kwargs)

        # Process crs
        if gdf.crs is None:
            warnings.warn('Boundaries object assuming external gdf data is in {} because no crs is given'.format(self.data.crs))
        elif gdf.crs != self.data.crs and self.data.crs is not None:
            warnings.warn('Boundaries object forced to transform external gdf data to new coordinate system')
            gdf = gdf.to_crs(self.data.crs)
        return gdf


    def fast_intersection(self, gdf, geometry_column = 'geometry', **kwargs):
        """

        Given a gdf of polygons or points, will calculate which boundary each polygon/point lies inside. Assumes that
        if the gdf is full of polygon data, each polygon will only intersect a single zip code.

        :param gdf: A geodataframe, in the same crs as the boundaries data.
        :param geometry_column: The geometry column of the gdf.
        :param kwargs: kwargs to pass to the underlying fast_polygon_intersection or points_intersect_multiple_polygons
            functions.
        :return: A pandas series mapping the index of the gdf to the indexes/names of the boundaries. If an index does not
            appear in the returned series, that is because the point/polygon for that index did not lie inside any of
            the boundaries polygons.

        """

        gdf = self.process_external_gdf(gdf, geometry_column = geometry_column)
        geometry_type = self.type_of_external_gdf(gdf)

        if geometry_type == 'Polygon':
            result = spatial_joins.fast_polygon_intersection(small_polygon_gdf = gdf, large_polygon_gdf = self.data,
                                                             small_geometry_column = geometry_column, **kwargs)
            return result
        elif geometry_type == 'Point':
            result = spatial_joins.points_intersect_multiple_polygons(points_gdf = gdf, polygons_gdf = self.data,
                                                                      points_geometry_column = geometry_column, **kwargs)
            return result
        else:
            raise TypeError('gdf geometry must be a shapely point or polygon')

    def push_features(self, gdf, features, geometry_column = 'geometry', **kwargs):

        """

        Given features in boundaries data, this "pushes" or calculates those features for another gdf of other (usually
        smaller) polygons. This is effectively a slower alternative to the fast_intersection method.

        :param GeoDataFrame gdf: A geodataframe of polygon geometry.
        :param features: The factors to calculate
        :param geometry_column: The geometry column of the gdf.
        :param kwargs: kwargs to pass to the underlying get_averages_by_area function.
        :return: The gdf but with 'factors' columns.

        """

        gdf = self.process_external_gdf(gdf, geometry_column = geometry_column)

        result = spatial_joins.get_averages_by_area(data_source = self.data, other_geometries = gdf, features = features,
                                                    other_geometries_column = geometry_column, **kwargs)
        return result

    def pull_features(self, gdf, features, geometry_column = 'geometry', **kwargs):

        """

        Given features in a non-boundaries gdf, this "pulls" or calculates those features for the boundary data attribute.
        This a slower alternative to the fast_intersection method.

        :param GeoDataFrame gdf: A geodataframe of polygon geometry.
        :param features: The factors to calculate from the gdf.
        :param geometry_column: The geometry column of the gdf.
        :param kwargs: kwargs to pass to the underlying get_averages_by_area function.
        :return: None but modifies self.data so that it has 'factors' columns.

        """

        container = self.data

        self.data = spatial_joins.get_averages_by_area(data_source = gdf, other_geometries = self.data, features = features,
                                                    other_geometries_column = geometry_column, **kwargs)

    def plot(self, **kwargs):
        self.data.plot(**kwargs)

    def __str__(self, **kwargs):
        return self.data.__str__(**kwargs)


class ZipBoundaries(Boundaries):
    """ Class for Zipcode Boundaries Data, wraps Boundaries class.

     :param ziplist: A list of zipcodes to subset to."""

    def __init__(self, ziplist = None, bounding_counties = None, bounding_polygon = None):

        # Put zips as strings
        if ziplist is not None:
            ziplist = list(map(str, ziplist))
            subset_by = 'index'
        else:
            subset_by = None
        self.ziplist = ziplist

        if isinstance(bounding_counties, str):
            bounding_counties = [bounding_counties]


        super().__init__(path = zip_boundaries_path, index_col = 'ZCTA5CE10', index_type = 'str', subset_by = subset_by,
                         subset_to = self.ziplist, bounding_counties=bounding_counties, bounding_polygon=bounding_polygon)

    def add_property_data(self, property_input, features = None, **kwargs):
        """ Add demand features from realtor/zillow inputs

        :param demand_input: A property_input class
        :param features: Optional features to subset to
        :param kwargs: **kwargs to pass to the .join method. Often useful to pass a rsuffix to prevent ValueErrors.
        :return: None, but changes self.data
        """

        property_data, metadata = property_input.process_property_data(features = features, geo_filter_values = 'all')
        property_data.index = list(property_data.index.map(str))
        self.data = self.data.join(property_data, how = 'left', **kwargs)


class BlockBoundaries(Boundaries):
    """ Class for block data, wraps Boundaries class, with a substantially different init method. After initialization,
    self.data is a gdf with a 'geometry' column, a variety of block data columns, and if cities is not None, it has a
    'place' column as well.

    :param data_layers: Iterable of codes for the data layer of the geodatabase.
    :param cities: City name (i.e. 'Austin'), or iterable of city names to filter by (i.e. ['Austin', 'Dallas'].
        Need to be in Texas or this won't work. Defaults to None.
   """

    def __init__(self, data_layers, cities = None, get_percent_residential = True):

        # Merge layers from geodatabase - this is a pretty cheap operation, relatively speaking (10 sec or so)
        geodata = gpd.read_file(texas_blocks_path, layer='ACS_2016_5YR_BG_48_TEXAS')
        geodata['geometry'] = geodata['geometry'].apply(lambda x: x[0])
        geodata['real_geometry_column'] = geodata['geometry'].copy()
        for data_layer in data_layers:
            data = gpd.read_file(texas_blocks_path, layer=data_layer)
            geodata = geodata.merge(data, how='inner', left_on='GEOID_Data', right_on='GEOID')
            geodata.rename({'GEOID_y': 'GEOID'}, axis='columns', inplace=True)

        geodata.set_index('GEOID', inplace=True)

        # Get the percent of land which is zoned residential inside the city limits
        if get_percent_residential:
            if os.path.exists(bg_percent_residential_path):
                percent_residential = pd.read_csv(bg_percent_residential_path, index_col=0).fillna(1)
                percent_residential.columns = ['percent_residential']
                geodata = geodata.join(percent_residential)
            else:
                warnings.warn('Cannot retrieve percent_residential for block data because percent_residential data has not been cached')

        geodata = gpd.GeoDataFrame(data=geodata[[x for x in geodata.columns if x != 'real_geometry_column']],
                                   geometry=geodata['real_geometry_column']) # This automatically renames the columns
        geodata.crs = {'init': 'epsg:4326'}

        if cities is not None:

            if isinstance(cities, str):
                cities = [cities]

            # Only consider blocks in cities
            place_shapes = gpd.read_file(texas_places_path)
            place_shapes = place_shapes.loc[place_shapes['NAME'].isin(cities)]
            mapper = spatial_joins.fast_polygon_intersection(geodata, place_shapes, large_name_column = 'NAME')
            geodata['place'] = geodata.index.map(mapper)
            geodata = geodata.loc[geodata['place'].notnull()]

        self.data = geodata
        self.sindex = geodata.sindex

def calc_percent_residential_in_block_groups():
    """ Calculates the percent zoned as residential in each block group. Should be run in the setup.py. This has been
    tested for accuracy but is not in the unittests.
    """

    from . import zoning

    # Get block and zoning data
    block_data = BlockBoundaries(data_layers = ['X19_INCOME'], cities = ['Austin', 'Dallas'], get_percent_residential = False)
    austin_block_data = block_data.data.loc[block_data.data['place'] == 'Austin']
    dallas_block_data = block_data.data.loc[block_data.data['place'] == 'Dallas']

    # Get zone data and simplify for performance. Note zone data is a geoseries not a geodataframe.
    dallas_zones = zoning.dallas_inputs.process_zoning_shapefile()
    dallas_zones = simple.process_geometry(dallas_zones.loc[dallas_zones['broad_zone'].isin(['Single Family', 'Multifamily'])])
    austin_zones = zoning.austin_inputs.process_zoning_shapefile()
    austin_zones = simple.process_geometry(austin_zones.loc[austin_zones['broad_zone'].isin(['Single Family', 'Multifamily'])])

    # Now find intersections
    for cityname, block_data, zone_data in zip(['Dallas', 'Austin'], [dallas_block_data, austin_block_data], [dallas_zones, austin_zones]):
        print('Starting to work on {}'.format(cityname))
        spatial_index = zone_data.sindex
        block_data['area'] = block_data['geometry'].area

        def calc_percent_res(polygon):
            percent_residential = spatial_joins.polygons_intersect_single_polygon(zone_data, polygon, spatial_index, factors = None, account_for_area = True, horiz = 1, vert = 1)
            return percent_residential

        block_data['percent_residential'] = block_data['geometry'].apply(calc_percent_res)

    results = pd.concat([austin_block_data['percent_residential'], dallas_block_data['percent_residential']])
    results.to_csv(bg_percent_residential_path)

austin_zips = [73301, 73344, 78704, 78705, 78708, 78713, 78714, 78715, 78701, 78702, 78703, 78709, 78710, 78711, 78712,
78716, 78717, 78718, 78719, 78720, 78721, 78722, 78723, 78728, 78729, 78730, 78731, 78734, 78724, 78725, 78726, 78727,
78732, 78733, 78735, 78736, 78739, 78741, 78745, 78746, 78749, 78750, 78751, 78752, 78755, 78756, 78757, 78758, 78759,
78737, 78738, 78742, 78744, 78747, 78748, 78753, 78754, 78760, 78761, 78762, 78763, 78769, 78772, 78773, 78780, 78781,
78799, 78764, 78765, 78766, 78767, 78768, 78774, 78778, 78779, 78783, 78785, 78789]
austin_zips = [str(z) for z in austin_zips]
dallas_zips = [75203, 75204, 75205, 75208, 75209, 75210, 75211, 75212, 75214, 75201, 75202, 75206, 75207, 75215, 75216,
75217, 75218, 75222, 75223, 75224, 75225, 75230, 75231, 75232, 75233, 75236, 75237, 75238, 75240, 75246, 75251, 75252,
75253, 75219, 75220, 75221, 75226, 75227, 75228, 75229, 75234, 75235, 75241, 75242, 75243, 75244, 75247, 75248, 75249,
75250, 75254, 75260, 75261, 75265, 75266, 75267, 75270, 75275, 75284, 75285, 75262, 75263, 75264, 75277, 75283, 75287,
75301, 75315, 75320, 75326, 75355, 75356, 75357, 75358, 75370, 75371, 75372, 75373, 75381, 75382, 75389, 75390, 75393,
75394, 75395, 75303, 75312, 75313, 75336, 75339, 75342, 75354, 75359, 75360, 75367, 75368, 75374, 75376, 75378, 75379,
75380, 75391, 75392, 75397, 75398]
dallas_zips = [str(z) for z in dallas_zips]
houston_zips = [77002, 77003, 77004, 77005, 77006, 77007, 77008, 77009, 77010, 77011, 77012, 77013, 77014, 77015, 77016,
                77017, 77018, 77019, 77020, 77021, 77022, 77023, 77024, 77025, 77026, 77027, 77028, 77029, 77030, 77031,
                77032, 77033, 77034, 77035, 77036, 77037, 77038, 77039, 77040, 77041, 77042, 77043, 77044, 77045, 77046,
                77047, 77048, 77049, 77050, 77051, 77053, 77054, 77055, 77056, 77057, 77058, 77059, 77060, 77061, 77062,
                77063, 77064, 77065, 77066, 77067, 77068, 77069, 77070, 77071, 77072, 77073, 77074, 77075, 77076, 77077,
                77078, 77079, 77080, 77081, 77082, 77083, 77084, 77085, 77086, 77087, 77088, 77089, 77090, 77091, 77092,
                77093, 77094, 77095, 77096, 77098, 77099, 77201, 77336, 77338, 77339, 77345, 77346, 77357, 77365, 77373,
                77375, 77377, 77379, 77386, 77388, 77396, 77401, 77406, 77407, 77429, 77433, 77447, 77449, 77450, 77477,
                77478, 77484, 77489, 77493, 77494, 77498, 77503, 77504, 77506, 77520, 77530, 77532, 77536, 77546, 77547,
                77571, 77587, 77598]
houston_zips = [str(z) for z in houston_zips]
