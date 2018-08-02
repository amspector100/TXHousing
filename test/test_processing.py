"""Tests the data_processing package."""

import os
from TXHousing.data_processing import zoning, property, parcel, boundaries
import TXHousing.chdir # This changes the directory to the parent directory
import matplotlib.pyplot as plt

import unittest

comprehensive_flag = False
while comprehensive_flag not in ['y', 'n']:
    comprehensive_flag = input("""For test_processing module, would you like to run a comprehensive test? Enter one of [Y/N]""").lower()


class TestZoningProcessing(unittest.TestCase):
    def setUp(self):
        print('Testing zoning module')

    # Test existence of data
    def test_data_existence(self):
        self.assertTrue(os.path.exists(zoning.north_texas_inputs.path),
                        'North Texas Zoning is not in the data directory; download it from http://data-nctcoggis.opendata.arcgis.com/datasets/2015-land-use')
        self.assertTrue(os.path.exists(zoning.austin_inputs.path),
                        'Austin zoning is not in the data directory; download it from https://data.austintexas.gov/Locations-and-Maps/Zoning/5rzy-nm5e')
        self.assertTrue(os.path.exists(zoning.dallas_inputs.path),
                        'Dallas zoning is not in the data directory; download it from https://dallasopendata.com/dataset/FY-2017-Base-Zoning/2wwb-qgic')
        self.assertTrue(os.path.exists(zoning.austin_regulations_path),
                        'Austin regulations data is not in the shared_data directory. Pull it from https://github.com/amspector100/TXHousing/tree/master/shared_data')
        self.assertTrue(os.path.exists(zoning.dallas_regulations_path),
                        'Dallas regulations data is not in the shared_data directory. Pull it from https://github.com/amspector100/TXHousing/tree/master/shared_data')

    # Test processing function
    def test_processing(self):
        try:
            dallas_zoning = zoning.dallas_inputs.process_zoning_shapefile(quietly = True, regulation_features = ['far', 'height', 'coverage'])
        except Exception as e:
            self.fail("process_zoning_shapefile unexpectedly raised {}".format(e))
        self.assertEqual(dallas_zoning.crs, {'init':'epsg:4326'}, 'Failed to transform or process crs correctly')

class TestPropertyProcessing(unittest.TestCase):
    def setUp(self):
        print('Testing property module')

    # Test existence of absolutely necessary data.
    def test_data_existence(self):
        self.assertTrue(os.path.exists(property.realtor_hotness_data.path),
                        'Realtor Hotness Data is not in the data directory; download it from https://realtor.com/research/data/')
        self.assertTrue(os.path.exists(property.realtor_core_inventory_sf.path),
                        'Realtor Core SF Inventory is not in the data directory; download it from https://realtor.com/research/data/')
        self.assertTrue(os.path.exists(property.realtor_core_inventory_mf.path),
                        'Realtor Core MF Inventory is not in the data directory; download it from https://realtor.com/research/data/')

    # Test processing function
    def test_processing(self):
        try:
            data, metadata = property.realtor_core_inventory.process_property_data()
        except Exception as e:
            self.fail("process_property_data unexpectedly raised {} for Realtor data".format(e))
        try:
            data, metadata = property.sfhomes_nbhd.process_property_data()
        except Exception as e:
            self.fail('process_property_data unexpectedly raised {} for Zillow data'.format(e))

# Note that this only tests ParcelProcessing for core municipalities
class TestParcelProcessing(unittest.TestCase):
    def setUp(self):
        print('Testing parcel module')

    # Test existence of necessary data
    def test_data_existence(self):
        self.assertTrue(os.path.exists(parcel.austin_parcel_path),
                        'Austin municipal parcel data is not in the data directory; download it from https://data.austintexas.gov/Building-and-Development/Land-Database-2016/nuca-fzpt')
        self.assertTrue(os.path.exists(parcel.dallas_parcel_data_path_2016),
                       'Dallas municipal parcel data is not in the data directory; download it from https://dallasopendata.com/Geography-Boundaries/Parcel-Shapefile/hy5f-5hrv')
        self.assertTrue(os.path.exists(parcel.harris_parcel_path_2018),
                        'Houston municipal parcel shapefile is not in the data directory; download it from http://pdata.hcad.org/download/index.html')
        self.assertTrue(os.path.exists(parcel.harris_parcel_building_res_path_2018),
                        'Houston municipal parcel data to merge with the shapefile is not in the data directory; download it from http://pdata.hcad.org/download/index.html')

    # Currently this takes several minutes - maybe I'll figure out a good way to do this in the future.
    def test_processing(self):
        if comprehensive_flag == 'y':
            try:
                houston_parcel_data = parcel.process_houston_parcel_data()
                del houston_parcel_data #( Save memory )
            except Exception as e:
                self.fail('process_houston_parcel_data unexpectedly raised {}'.format(e))

class TestBoundariesProcessing(unittest.TestCase):

    def setUp(self):
        print('Testing boundaries module')

    def test_data_existence(self):
        self.assertTrue(os.path.exists(boundaries.zip_boundaries_path),
                        'Zipcode boundaries data is not in the data directory; download it from https://census.gov/geo/maps-data/data/cbf/cbf_zcta.html')
        self.assertTrue(os.path.exists(boundaries.county_boundaries_path),
                        'County boundaries data is not in the data directory; download it from https://census.gov/geo/maps-data/data/cbf/cbf_counties.html')
        self.assertTrue(os.path.exists(boundaries.texas_blocks_path),
                        'Texas block data is not in the data directory; download it from https://census.gov/geo/maps-data/data/tiger-data.html; use 2012-2016 selected tables')



    def test_boundaries_parent_class(self):

        # This actually has remarkably high code coverage for __init__ although it might not seem like it initially
        try:
            zipdata = boundaries.ZipBoundaries(ziplist = None, bounding_counties = ['Harris', 'Travis'])
        except Exception as e:
            self.fail('Boundaries class initialization unexpectedly raised {}'.format(e))

        try:
            counties = boundaries.Boundaries(path = boundaries.county_boundaries_path,
                                             subset_by = 'STATEFP',
                                             subset_to = ['48'])
            result = zipdata.fast_intersection(counties.data)
            print(result)
        except Exception as e:
            self.fail('Boundaries.fast_intersection unexpectedly raised {}'.format(e))



if __name__ == '__main__':
    unittest.main()