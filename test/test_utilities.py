"""Tests core utilities. This is NOT a comprehensive testing suite."""

from TXHousing.utilities import simple, spatial_joins, measurements
import unittest
import warnings

# For testing
import numpy as np
import shapely.geometry
import shapely.affinity
import pandas as pd
import geopandas as gpd

# Set random state for repeatability
np.random.seed(110)

class TestSimpleUtilities(unittest.TestCase):

    def setUp(self):

        # Create points/features
        x_coords = np.random.randn(100)
        y_coords = np.random.randn(100)
        self.features = np.random.randn(100)
        self.points = [shapely.geometry.point.Point(x, y) for x, y in zip(x_coords, y_coords)]

        # Create polygons
        self.polygon = shapely.geometry.polygon.Polygon([[0, 1], [0.5, 0.5], [1.5, 0.5], [0.5, -0.5], [3, -0.5], [1, 3]])

    def test_fragment(self):
        try:
            grid = simple.fragment(self.polygon)
        except:
            print('hello')
        self.assertAlmostEqual(self.polygon.area, sum(p.area for p in grid), places = 6, msg = 'simple.fragment failed to precisely partition polygon area')
        self.assertTrue(np.all(p.area < self.polygon.area/100 for p in grid), 'simple.fragment failed to sufficiently decrease the size of fragments')

    def test_point_grid(self):
        # Test that point_grid works and that its precision is decent
        gdf = gpd.GeoDataFrame(data = pd.DataFrame(pd.Series(self.features), columns = ['value']), geometry = self.points)
        try:
            result = simple.make_point_grid(gdf, factor = 'value', by = 'mean')
            gdf_area = shapely.geometry.multipoint.MultiPoint(gdf['geometry'].values).convex_hull.area
            self.assertAlmostEqual(result['geometry'].area.sum(), gdf_area, places = 6, msg = 'simple.make_point_grid insufficiently precise in partitioning areas')
        except Exception as e:
            self.fail('simple.make_point_grid unexpectedly raised {} for by = "mean"'.format(e))
        try:
            result = simple.make_point_grid(gdf, factor = None, by = 'mean')
        except Exception as e:
            self.fail('simple.make_point_grid unexpectedly raised {} for factors = None'.format(e))

    def test_urban_core(self):

        # Check to make sure urban core is the right area
        lat = 30.2331
        long = -97.3421
        radius = 100
        urban_core = simple.get_urban_core(lat, long, radius = radius, scale = 1)
        urban_core.crs = {'init':'epsg:4326'}
        urban_core = urban_core.to_crs({'init':'epsg:2276'})
        self.assertAlmostEqual(urban_core.area.sum(), np.pi * radius**2, places = -2, msg = 'simple.get_urban_core insufficiently precise')

    def test_hex_color(self):
        color = np.array([240,248,255])/255
        self.assertEqual(simple.convert_to_hex(color), '#F0F8Ff', 'simple.convert_to_hex fails to properly convert colors')


class TestSpatialJoins(unittest.TestCase):

    def setUp(self):

        # Generate geometries
        self.p1 = shapely.geometry.polygon.Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])
        self.p2 = shapely.geometry.polygon.Polygon([[2, 0], [1.5, -1], [1, 0]])
        self.p3 = shapely.geometry.polygon.Polygon([[0, 0], [0, 0.6], [1, 1.6], [1, 0]])

        # Large polygons
        self.large_polygons = gpd.GeoDataFrame(data = pd.DataFrame(pd.Series(['name1', 'name2']), columns = ['name']), geometry = [self.p1, self.p2])


        # Create version of large_polygons that will throw a warning because a point will lie in multiple polygons;
        # then add noncategorical data to the initial large_polygons instance.
        self.large_polygons_v2 = self.large_polygons.copy()
        self.large_polygons_v2.loc[self.large_polygons_v2.index[-1] + 1] = ['name3', self.p3]

        self.large_polygons['value'] = np.arange(0, self.large_polygons.shape[0], 1) + 2
        self.large_polygons['value2'] = np.arange(0, self.large_polygons.shape[0], 1) + 3


        # Points
        self.points = [shapely.geometry.point.Point([0.5, 0.5]), shapely.geometry.point.Point([0.3, 0.4]),
                       shapely.geometry.point.Point(1.3, 0.3), shapely.geometry.point.Point([0.4, 0.2])]
        self.points_gdf = gpd.GeoDataFrame(data = pd.DataFrame(pd.Series([1,2,3,4]), columns = ['value']), geometry = self.points)
        self.points_gdf['value2'] = [0, 1, 2, 3]

        # Create small gdf
        small_gdf = self.points_gdf.copy()
        small_gdf['geometry'] = small_gdf['geometry'].apply(lambda x: x.buffer(5))
        self.small_gdf = small_gdf


    def test_points_intersections(self):

        # Check that the intersections are correct for points_intersect_multiple_polygons
        result = spatial_joins.points_intersect_multiple_polygons(points_gdf = self.points_gdf,
                                                                  polygons_gdf = self.large_polygons)
        self.assertEqual(set(result.index.tolist()), set([1, 3, 0]),
                         'spatial_joins.points_intersect_multiple_polygons does not correctly check intersections')
        self.assertEqual(result.unique().tolist(), [0],
                         'spatial_joins.points_intersect_polygons does not correctly check intersections')

        # Check that calculations and intersections are correct for points_intersect_single_polygon
        categorical_result = spatial_joins.points_intersect_single_polygon(points = self.points_gdf, polygon = self.p1,
                                                                           spatial_index = self.points_gdf.sindex,
                                                                           factors = ['value', 'value2'], categorical = True)

        self.assertEqual(np.unique(categorical_result.values).tolist(), [1],
                         'spatial_joins.points_intersect_single_polygon incorrectly calculates intersections')

        self.assertEqual(categorical_result.index.tolist(), [(1, 0), (2, 1), (4, 3)],
                         """spatial_joins.points_intersect_single_polygon either incorrectly calclulates intersections 
                         or incorrectly groups by categorical factors""")

        continuous_result = spatial_joins.points_intersect_single_polygon(points = self.points_gdf, polygon = self.p1,
                                                                           spatial_index = self.points_gdf.sindex,
                                                                           factors = ['value', 'value2'], categorical = False, by = 'median')

        self.assertIsInstance(continuous_result, pd.Series,
                              'spatial_joins.points_intersect_single_polygon does not return a series for multiple non-categorical factors')
        self.assertEqual(continuous_result['value'], 2,
                         """spatial_joins.points_intersect_single_polygon  either incorrectly calclulates intersections 
                         or incorrectly calculates medians""")

    def test_polygons_intersection(self):

        # Check to make sure that fast_polygon_intersection and points_intersect_multiple_polygons properly warn the user
        # when points/small polygons intersect multiple large polygons.

        with warnings.catch_warnings(record = True) as w:

            # Ensures all warnings are triggered
            warnings.simplefilter("always")

            # Run fast_polygon intersection and see if we get the right intersections although the names might be different
            result = spatial_joins.fast_polygon_intersection(small_polygon_gdf = self.small_gdf,
                                                             large_polygon_gdf = self.large_polygons_v2,
                                                             large_name_column = 'name')
            self.assertEqual(set(result.index.tolist()), set([0,1,3]),
                             'spatial_joins.fast_polygon_intersection incorrectly calculates intersections')


            assert len(w) == 1
            assert issubclass(w[-1].category, UserWarning)
            assert 'only returns one polygon per point' in str(w[-1].message)


        # Now test polygons_intersect_single_polygon
        spatial_index = self.large_polygons.sindex
        categorical_result = spatial_joins.polygons_intersect_single_polygon(small_polygons = self.large_polygons,
                                                                             polygon = self.p3,
                                                                             spatial_index = spatial_index,
                                                                             factors = 'name',
                                                                             categorical = True)
        self.assertAlmostEqual(categorical_result['name2'], 0,
                               msg = 'spatial_joins.polygons_intersect_single_polygon incorrectly calculates intersections')

        continuous_result = spatial_joins.polygons_intersect_single_polygon(small_polygons = self.large_polygons,
                                                                            polygon = self.p3,
                                                                            spatial_index = spatial_index,
                                                                            factors = ['value', 'value2'],
                                                                            by = 'mean',
                                                                            categorical = False)

        expected_result = 2*(self.large_polygons.loc[0, 'geometry'].intersection(self.p3).area/self.p3.area)
        self.assertAlmostEqual(continuous_result['value'], expected_result,
                               msg = 'spatial_joins.polygons_intersect_single_polygon incorrectly calculates area-weighted means')




class TestMeasurements(unittest.TestCase):

    def setUp(self):
        self.HK = shapely.geometry.point.Point([114.1095, 22.3964])  # Hong Kong coordinates
        self.Pelham = shapely.geometry.point.Point([-73.8079, 40.9098])  # Pelham, NY coordinates
        self.NYC = shapely.geometry.point.Point([-74.0060, 40.7128])  # NYC coordinates

        # Create some polygons for testing
        p1 = shapely.geometry.polygon.Polygon([[0, 1], [0.5, 0.5], [1.5, 0.5], [0.5, -0.5], [3, -0.5], [1, 3]])
        p2 = shapely.geometry.polygon.Polygon([[2, 0], [1.5, -1], [1, 0]])
        p3 = shapely.geometry.polygon.Polygon([[10, 0], [10.5, -1], [9, -3]])
        p4 = shapely.geometry.polygon.Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])
        self.gdf = gpd.GeoDataFrame(data = pd.DataFrame(pd.Series([0, 4, 2, 2]), columns = ['values']), geometry = [p1, p2, p3, p4])
        self.gdf.crs = {'init':'epsg:2277'}

    def test_data_existence(self):

        texas_places_path


    def test_haversine(self):
        NY_to_HK = measurements.haversine(point1=self.HK, point2=self.NYC)
        self.assertTrue(abs(NY_to_HK - 8046) < 10,
                        'measurements.haversine function is insufficiently precise for long distances')
        NY_to_Pelham = measurements.haversine(point1=self.NYC, point2=self.Pelham)
        self.assertTrue(abs(NY_to_Pelham - 17.106) < 0.03,
                        'measurements.haversine function is insufficiently precise for short distances')

    def test_order_radii(self):

        data = pd.DataFrame(np.arange(12).reshape(6, 2))
        data.index = [1, '2', '3', '4', '5', '10+']
        data = measurements.order_radii(data)
        self.assertTrue(data.index[0] < data.index[-1], 'measurements.order_radii fails to properly order radii')

    def test_point_rings(self):

        gdf = self.gdf.copy()
        gdf['geometry'] = gdf['geometry'].apply(lambda x: shapely.affinity.scale(x, xfact = 10000, yfact = 10000))
        gdf = gdf.to_crs({'init':'epsg:4326'})
        gdf['centroids'] = gdf['geometry'].centroid
        gdf = gdf.set_geometry('centroids')

        try:
            result = measurements.points_intersect_rings(gdf, lat=3.48, long=-105.9, factor='values', step=0.5,
                                                         categorical=True, by='mean',
                                                         geometry_column='centroids', per_square_mile=False, maximum=6.5)
        except Exception as e:
            self.fail('measurements.points_intersect_rings unexpectedly raised {}'.format(e))
        self.assertEqual(result.iloc[0].sum(), 3)


    def test_polygon_rings(self):

        self.gdf = self.gdf.to_crs({'init':'epsg:4326'})
        self.lat_center = 3.442338
        self.long_center = -105.98319

        # Start by testing polygons intersect rings for categorical data
        try:
            result1 = measurements.polygons_intersect_rings(self.gdf, lat = self.lat_center, long = self.long_center, factor = 'values')
        except Exception as e:
            self.fail('measurements.polygons_intersect_rings unexpectedly raised {}'.format(e))

        # Handle other tests in case we could not actually generate result1
        try:
            self.assertAlmostEqual(result1[0].notnull().sum(), 1, places = 3,
                                   msg = 'measurements.polygons_intersect_rings incorrectly calculates distance from center')
            self.assertAlmostEqual(result1.iloc[0].sum(), 1, places = 3,
                                   msg = 'measurements.polygons_intersect_rings does not properly calculate area percentages')
        except NameError:
            pass

        # Now polygons intersect rings for continuous data. Just make sure it works.
        try:
            result2 = measurements.polygons_intersect_rings(self.gdf, lat = self.lat_center, long = self.long_center, factor = 'values', step = 0.001, maximum = 0.003, categorical = False)
        except Exception as e:
            self.fail('measurements.polygons_intersect_rings unexpectedly raised {} for noncategorical data'.format(e))

        # Handle other tests in case we could not generate result1 or 2
        try:
            expected_value = 2*result1.iloc[0, 2] + 4*result1.iloc[0, 1]
            self.assertAlmostEqual(result2.iloc[0], expected_value, places = 1,
                                   msg = 'measurements.polygons_intersect_rings does not properly calculate area-weighted means')
        except NameError:
            pass

    def test_area_calculations(self):
        gdf = measurements.get_area_in_units(self.gdf, scale = 1)
        self.assertAlmostEqual(gdf['area'].iloc[3], 1, places = 5, msg = 'measurements.get_area_in_units is not precise enough')

if __name__ == '__main__':
    unittest.main()
