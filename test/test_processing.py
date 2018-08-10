"""Tests the data_processing package."""

import os
from TXHousing.data_processing import zoning, property, parcel, permit, boundaries
import TXHousing.chdir # This changes the directory to the parent directory

import unittest
import warnings
import numpy as np

# Set random state for repeatability
np.random.seed(110)

comprehensive_flag = False
while comprehensive_flag not in ['y', 'n']:
    comprehensive_flag = input("""For test_processing module, would you like to run a comprehensive test? Enter one of [Y/N]""").lower()

class TestZoningProcessing(unittest.TestCase):
    def setUp(self):
        pass

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
        self.assertTrue(os.path.exists(zoning.houston_spec_min_lot_path),
                        'Houston special minimum lot zone shapefile is not in the data directory. Download it from https://cohgis-mycity.opendata.arcgis.com/datasets/minimum-lot-size')
        self.assertTrue(os.path.exists(zoning.houston_spec_min_lot_path),
                        'Houston special setback zone shapefile is not in the data directory. Download it from https://cohgis-mycity.opendata.arcgis.com/datasets/special-minimum-building-lines')

    def test_historic_district_data_existence(self):
        self.assertTrue(os.path.exists(zoning.tx_hd_path),
                        'Texas national historic landmarks data is not in the data directory. Download it from https://atlas.thc.state.tx.us/Data/GISData')
        self.assertTrue(os.path.exists(zoning.austin_landmark_path),
                        'Austin historic landmarks data is not in the data directory. Download it from https://data.austintexas.gov/Locations-and-Maps/Historical-Landmarks/vvuz-m3y4')
        self.assertTrue(os.path.exists(zoning.dallas_historic_overlay_path),
                        'Dallas historic overlay data is not in the data directory. Download it from https://gis.dallascityhall.com/shapefileDownload.aspx')
        self.assertTrue(os.path.exists(zoning.dallas_historic_subdistricts_path),
                        'Dallas historic subdistricts data is not in the data directory. Download it from https://gis.dallascityhall.com/shapefileDownload.aspx')
        self.assertTrue(os.path.exists(zoning.houston_historic_districts_path),
                        'Houston historic districts data is not in the data directory. It is currently unavailable online. Email the package maintainer to get it.')
        self.assertTrue(os.path.exists(zoning.houston_historic_landmarks_path),
                        'Houston historic landmarks data is not in the data directory. It is currently unavailable online. Email the package maintainer to get it.')


    # Test processing function
    def test_processing(self):
        try:
            dallas_zoning = zoning.dallas_inputs.process_zoning_shapefile(quietly = True, regulation_features = ['far', 'height', 'coverage'])
        except Exception as e:
            self.fail("process_zoning_shapefile unexpectedly raised {}".format(e))
        self.assertEqual(dallas_zoning.crs, {'init':'epsg:4326'}, 'Failed to transform or process crs correctly')

class TestPropertyProcessing(unittest.TestCase):
    def setUp(self):
        pass

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
            self.fail("process_property_data, for Realtor data, unexpectedly raised {}".format(e))
        try:
            data, metadata = property.sfhomes_nbhd.process_property_data()
        except Exception as e:
            self.fail('process_property_data, for Zillow data, unexpectedly raised {}'.format(e))

# Note that this only tests ParcelProcessing for core municipalities
class TestParcelProcessing(unittest.TestCase):
    def setUp(self):
        pass

    # Test existence of necessary data
    def test_data_existence(self):
        self.assertTrue(os.path.exists(parcel.austin_parcel_path),
                        'Austin municipal parcel data is not in the data directory; download it from https://data.austintexas.gov/Building-and-Development/Land-Database-2016/nuca-fzpt')
        self.assertTrue(os.path.exists(parcel.dallas_county_parcel_path),
                        """Dallas county parcel data (used to calc municipal parcel data) is not in the data directory; 
                        download it from http://dallascad.org/contact.aspx""")
        self.assertTrue((os.path.exists(parcel.dallas_county_appraisal_path) & os.path.exists(parcel.dallas_county_res_path)),
                        """Dallas county parcel data features are not in the data directory download them from
                         http://dallascad.org/contact.aspx""")
        self.assertTrue(os.path.exists(parcel.harris_parcel_path_2018),
                        'Houston municipal parcel shapefile is not in the data directory; download it from http://pdata.hcad.org/download/index.html')
        self.assertTrue(os.path.exists(parcel.harris_parcel_building_res_path_2018),
                        'Houston municipal parcel data to merge with the shapefile is not in the data directory; download it from http://pdata.hcad.org/download/index.html')

    # Currently this takes several minutes - maybe I'll figure out a good way to do this in the future.
    def test_processing(self):
        if comprehensive_flag == 'y':
            try:
                houston_parcel_data = parcel.process_houston_parcel_data(county_level = False)
                del houston_parcel_data #( Save memory )
            except Exception as e:
                self.fail('process_houston_parcel_data unexpectedly raised {}'.format(e))

class TestBoundariesProcessing(unittest.TestCase):

    def setUp(self):
        pass

    def test_data_existence(self):
        self.assertTrue(os.path.exists(boundaries.zip_boundaries_path),
                        'Zipcode boundaries data is not in the data directory; download it from https://census.gov/geo/maps-data/data/cbf/cbf_zcta.html')
        self.assertTrue(os.path.exists(boundaries.county_boundaries_path),
                        'County boundaries data is not in the data directory; download it from https://census.gov/geo/maps-data/data/cbf/cbf_counties.html')
        self.assertTrue(os.path.exists(boundaries.texas_blocks_path),
                        'Texas block data is not in the data directory; download it from https://census.gov/geo/maps-data/data/tiger-data.html; use 2012-2016 selected tables')
        self.assertTrue(os.path.exists(boundaries.texas_places_path),
                        'Texas places boundary data is not in the data directory; download it from https://census.gov/geo/maps-data/data/cbf/cbf_place.html')
        self.assertTrue(os.path.exists(boundaries.ua_path),
                        'Texas urban areas boundary data is not in the data directory; download it from https://www.census.gov/geo/maps-data/data/cbf/cbf_ua.html')
        if os.path.exists(boundaries.bg_percent_residential_path) == False:
            warnings.warn("""Percent_residential for block data has not been cached. This may reduce the accuracy of a 
            couple of population calculations. Either run calculate_percent_residential in the boundaries module, or 
            simple pull the cache from https://github.com/amspector100/TXHousing/tree/reorganization/shared_data""")

    def test_boundaries_class(self):


        # This actually has remarkably high code coverage for __init__
        try:
            zipdata = boundaries.ZipBoundaries(ziplist = None, bounding_counties = ['Travis'])
        except Exception as e:
            self.fail('Boundaries class initialization unexpectedly raised {}'.format(e))

        if comprehensive_flag == 'y':
            try:
                counties = boundaries.Boundaries(path = boundaries.county_boundaries_path,
                                                 subset_by = 'STATEFP',
                                                 subset_to = ['48'])
                result = counties.fast_intersection(zipdata.data)
            except Exception as e:
                self.fail('Boundaries.fast_intersection unexpectedly raised {}'.format(e))

            try:
                self.assertEqual(result.unique().tolist(), [227, 692, 697, 2240, 2273, 2294, 2504],
                                 'boundaries.Boundaries.fast_intersection does not calculate correct intersections')
            except NameError:
                pass

        # Now test block geodata function
        try:
            blockdata = boundaries.BlockBoundaries(data_layers = ['X01_AGE_AND_SEX', 'X08_COMMUTING', 'X19_INCOME'], cities = ['Austin'])
        except Exception as e:
            self.fail('BlockBoundaries class initialization unexpectedly raised {}'.format(e))

        # Test pull and pushing data
        if comprehensive_flag == 'y':

            zipdata.data = zipdata.data.sample(15)
            zipdata.pull_features(blockdata.data, features = 'B01001e1', account_method = 'percent_residential')
            zipdata.data = blockdata.push_features(zipdata.data, features=['B01002e1'], account_method='water')

class TestPermitProcessing(unittest.TestCase):

    def setUp(self):
        pass

    def test_raw_data_existence(self):

        self.assertTrue(os.path.exists(permit.austin_permit_path),
                        'Austin permit data is not in the data directory; download it from https://data.austintexas.gov/Building-and-Development/Issued-Construction-Permits/3syk-w9eu')
        self.assertTrue(os.path.exists(permit.dallas_permit_path),
                        'Dallas permit data is not in the data directory; download it from https://www.dallasopendata.com/City-Services/Building-Inspection-Master-Permits/ks9j-qkj8')
        self.assertTrue(os.path.exists(permit.houston_structural_permits_path),
                        'Houston structural permit data is not in the data directory; download it from https://cohgis-mycity.opendata.arcgis.com/datasets/permits-wm-structural')
        self.assertTrue(os.path.exists(permit.houston_structural_permits_path),
                        'Houston demolition permit data is not in the data directory; download it from https://cohgis-mycity.opendata.arcgis.com/datasets/demolition-ilms-code-sd-1')

    def test_processed_data_existence(self):

        self.assertTrue(os.path.exists(permit.dpm_save_path),
                        """Dallas corrected permit data is not in the data directory; pull it from https://github.com/amspector100/TXHousing/tree/reorganization/shared_data, 
                        or run permit.correct_dallas_permit_data()""")
        self.assertTrue(os.path.exists(permit.houston_permit_statuses_path),
                        """Houston permit updated approval statues are not in the data directory; pull them from https://github.com/amspector100/TXHousing/tree/reorganization/shared_data, 
                        or run permit.scrape_houston_permit_data()""")


    def test_basic_processing_functions(self):

        if comprehensive_flag == 'y':

            try:
                data = permit.process_austin_permit_data(searchfor = ['S.F.'])
            except Exception as e:
                self.fail('permit.process_austin_permit_data unexpectedly raised {}'.format(e))
            try:
                data = permit.get_corrected_dallas_permit_data()
            except Exception as e:
                self.fail('permit.get_corrected_dallas_permit_data unexpectedly raised {}'.format(e))
            try:
                data = permit.process_houston_permit_data()
            except Exception as e:
                self.fail('permit.get_corrected_dallas_permit_data unexpectedly raised {}'.format(e))

if __name__ == '__main__':
    unittest.main()