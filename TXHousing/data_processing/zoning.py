"""Processes zoning data for Austin, Dallas, and all of North Texas"""

import time
import pandas as pd
import geopandas as gpd
import warnings
from ..utilities import simple

class Zoning_Input():

    """ Simple class which binds information about zoning datasets. Also, the lat/long features define the city
        center for each city."""

    def __init__(self, path, feature, separator, base_zones, proj4string = None, lat = 0, long = 0, zoom = 10,
                 crs = None, regulations_path = None):
        self.path = path
        self.feature = feature
        self.separator = separator
        self.base_zones = base_zones
        self.proj4string = proj4string
        self.lat = lat
        self.long = long
        self.zoom = zoom
        self.crs = crs
        self.regulations_path = regulations_path

    def process_zoning_shapefile(self, overlay = {'-MU':'Multifamily'}, broaden = True, parse_base_zones = True, regulation_features = None,
                                 to_latlong = True, quietly = False):
        """
        :param broaden: Default true. If true, will decode zoning data into a broader classification (i.e. sf, mf)
            as specified by the zoning input class - this processed zoning data will go in a column labelled "broad_zone."
        :type broaden: Boolean
        :param overlay: A dictionary which maps overlay strings to broadened outputs. This is mostly useful for dealing with
            mixed use overlays in Austin (there are no mixed-use base zones, so passing in '-mu':'Multifamily'} will help
            the function recognize that mixed-use zones are multifamily zones.
        :param parse_base_zones: Boolean. Default true. If true, will use the regulations path provided by the input to
            parse the base zones of the zoning data.
        :param regulation_features: A list of features to retrieve from the regulation data. Default None.
        :param to_latlong: If True, will try to transform the data to lat/long.
        :param quietly: Default false. If true, will suppress all print statements and make the function "quiet"
        :type quietly: Boolean
        """

        # Begin timing
        time0 = time.time()

        # Read data and process zone codes
        if quietly is not True:
            print('Reading file')
        raw_data = gpd.read_file(self.path)
        raw_data = simple.process_geometry(raw_data, drop_multipolygons=False)
        if quietly is not True:
            print('Finished reading file, took {}'.format(time.time() - time0))

        if broaden:
            def get_zone(text):
                split_text = text.split(self.separator)[0]
                for key in self.base_zones:
                    if split_text in self.base_zones[key]:
                        return (key)
                return ([key for key in self.base_zones][-1])

            raw_data['broad_zone'] = raw_data[self.feature].apply(get_zone)

            # Now put account for overlays
            if overlay is not None:
                for key in [key for key in overlay]:
                    raw_data.loc[raw_data[self.feature].str.contains(key), 'broad_zone'] = overlay[key]

        # Parse base zones and retrieve regulatory data
        if parse_base_zones:

            if self.regulations_path is not None:

                # Start by parsing base zones ------------

                # Get regulations data - this has all the data and the list of base zones.
                reg_data = pd.read_csv(self.regulations_path, index_col=0, encoding='Latin-1')

                # If the data includes non-base zones, only include base zones in this particular search
                if 'class' in reg_data.columns:
                    reg_data = reg_data.loc[reg_data['class'] == 'base']

                reg_data.index = [str(ind) for ind in reg_data.index]

                # Process the zoning data - this function is well defined, I checked
                base_zones = reg_data.index.unique().tolist()

                def process_zone(text):
                    for i in base_zones:
                        if i in text:
                            return i
                    return 'Unknown'

                # Apply and get base zone
                raw_data['base_zone'] = raw_data[self.feature].apply(process_zone)

                # Now get regulation features ------------
                if regulation_features is not None:
                    if isinstance(regulation_features, str):
                        regulation_features = [regulation_features]
                    for regulation_feature in regulation_features:
                        raw_data[regulation_feature] = raw_data['base_zone'].map(reg_data[regulation_feature])

            else:
                warnings.warn(
                    """parse_base_zones = True, but a regulations path is required to parse base_zones from zoning data.""")

        # Manually set CRS if necessary, then transform
        if raw_data.crs is not None:
            raw_data.crs = self.crs

        if to_latlong:
            try:
                raw_data = raw_data.to_crs({'init': 'epsg:4326'})
            except ValueError as e:
                warnings.warn('Could not transform data to lat/long')

        if quietly is not True:
            print('Finished processing zones, mixed use overlay, and crs, took {}'.format(time.time() - time0))

        return raw_data

# Downloaded from http://data-nctcoggis.opendata.arcgis.com/datasets/2015-land-use
north_texas_inputs = Zoning_Input(path = 'data/Zoning/north_texas_land_use_2015/2015_Land_Use.shp',
                                   feature = 'CATEGORY',
                                   separator = 'no_separator',
                                   proj4string = 'EPSG:4326',
                                   base_zones = { 'Single Family':['Single family'],
                                               'Multifamily':['Mixed use', 'Multi-family'],
                                               'Other Residential':['Residential acreage'],
                                               'Other':[]},
                                   lat = 32.7767,
                                   long = -96.7970,
                                   zoom = 10,
                                  crs={'init': 'epsg:2276'})

# Downloaded from https://data.austintexas.gov/Locations-and-Maps/Zoning/5rzy-nm5e
austin_regulations_path = "shared_data/austin zoning standards.csv"
austin_inputs = Zoning_Input(path = "data/Zoning/austin_zoning/geo_export_571668ee-52f1-4ac9-a4e0-3b8bb348eae7.shp",
                              feature = 'zoning_zty',
                              separator = '-',
                              proj4string = 'EPSG:4326',
                              base_zones = {'Single Family':['SF'], 'Multifamily':['MF'],
                                            'Other Residential':['MH', 'RR', 'LA'], 'Other':[]},
                              lat = 30.267,
                              long = -97.743,
                              zoom = 9,
                              crs = {'init': 'epsg:4326'},
                              regulations_path=austin_regulations_path)

# From dallas's open data site
dallas_regulations_path = 'shared_data/dallas zoning standards.csv'
dallas_inputs = Zoning_Input(path = "data/Zoning/dallas_zoning/BaseZoning.shp",
                              feature = 'ZONE_DIST',
                              separator = '-',
                              proj4string = 'EPSG:2276',
                              base_zones = {'Single Family':['A(A)', 'D(A)', 'TH', 'R'],
                                            'Multifamily': ['CH', 'MF', 'MU'],
                                            'Other Residential':['MH(A)'],
                                            'Other':[]},
                              lat = north_texas_inputs.lat,
                              long = north_texas_inputs.long,
                              zoom = 9,
                              crs = {'init':'epsg:2276'},
                              regulations_path=dallas_regulations_path)

# Houston has no base zones but this is useful for making maps
houston_inputs = Zoning_Input(path = None,
                               feature = None,
                               separator = None,
                               proj4string = None,
                               base_zones = None,
                               lat = 29.7604,
                               long = -95.3698,
                               zoom = 9)

# Selected surrounding cities for Austin ---------------
# Some notes on methodology, for consistency. I classify:
# (a) Areas which allow both multifamily and single family dwellings as multifamily.
# (b) Areas which allow mixed-use development (i.e. commercial + residential) as multifamily in most cases.
# (c) If multifamily housing is permitted but with some restrictions, the area is classified as multifamily housing.


round_rock_inputs = Zoning_Input(path = "data/Zoning/round_rock_zoning/ZONING_1.shp",
                                 feature = 'BASE_ZONIN',
                                 separator = '-',
                                 base_zones = {'Single Family':['SF', 'TH', 'SF1', 'SF2', 'SF3', 'SF4', 'SF5', 'SF6'],
                                               'Multifamily': ['MF', 'MU', 'TF'],
                                               'Other Residential':[],
                                               'Other':[]},
                                 lat = austin_inputs.lat,
                                 long = austin_inputs.long,
                                 crs = {'init':'epsg:2277'})

pflugerville_inputs = Zoning_Input(path = "data/Zoning/pflugerville_zoning/Zoning_Districts.shp",
                                   feature = 'ZOINING_TY',
                                   separator = '-',
                                   base_zones = {'Single Family':['A', 'SF'],
                                                 'Multifamily': ['2', 'CL3', 'CL4', 'CL5'],
                                                 'Other Residential':['MH'],
                                                 'Other':[]},
                                   lat = austin_inputs.lat,
                                   long = austin_inputs.long,
                                   crs = {'init':'epsg:2277'})

georgetown_inputs = Zoning_Input(path = "data/Zoning/georgetown_zoning/Zoning.shp",
                                   feature = 'ZONE',
                                   separator = '-',
                                   base_zones = {'Single Family':['AG', 'RE', 'RL', 'RS', 'TH'],
                                                 'Multifamily': ['TF', 'MF', 'MU', 'MU-DT', 'MUDT'],
                                                 'Other Residential':['MH'],
                                                 'Other':[]},
                                   lat = austin_inputs.lat,
                                   long = austin_inputs.long,
                                   crs = {'init':'epsg:2277'})

cedar_park_inputs = Zoning_Input(path = "data/Zoning/cedar_park_zoning/Zoning__Zoning_Districts.shp",
                                   feature = 'ZoningType',
                                   separator = ' - ',
                                   base_zones = {'Single Family':['RA', 'SR', 'SU', 'UR'],
                                                 'Multifamily': ['MU', 'MF'],
                                                 'Other Residential':[],
                                                 'Other':[]},
                                   lat = austin_inputs.lat,
                                   long = austin_inputs.long,
                                   crs = {'init':'epsg:2277'})

hutto_inputs = Zoning_Input(path = "data/Zoning/hutto_zoning/Zoning_Districts.shp",
                            feature = 'ZONING',
                            separator = 'no_separator',
                            # Not totally sure about urban residential/residential classifications
                            base_zones = {'Single Family':['Single Family', 'Residential'],
                                          'Multifamily': ['Two Family', 'Multi-Family',
                                                          'Urban Residential', 'Co-op District'],
                                          'Other Residential':[],
                                          'Other':[]},
                            lat = austin_inputs.lat,
                            long = austin_inputs.long,
                            crs = {'init':'epsg:2277'})

def get_austin_surrounding_zones():
    """Based on the zoning inputs in this file, returns a geodataframe of zoning polygons for Austin AND select
     surrounding areas, with two columns: broad_zone and geometry. This is purely a convenience function which wraps a
     variety of Zoning_Input.process_zoning_shapefile calls. """

    data = gpd.GeoDataFrame(data = None, geometry = None)
    inputs = [round_rock_inputs, georgetown_inputs, pflugerville_inputs, hutto_inputs, cedar_park_inputs, austin_inputs]
    for input in inputs:
        to_add = input.process_zoning_shapefile()[['broad_zone', 'geometry']]
        data = pd.concat([data, to_add])
    return data

# Paths for historic districts - thankfully do not require much processing -----------------------------------------
tx_hd_path = "data/Zoning/national_register_hist_dists_shp/NationalRegisterPY.shp"
austin_landmark_path = "data/Zoning/austin_historical_landmarks/geo_export_ce453b58-d8ca-47d4-8934-76d636d24ca3.shp"

dallas_historic_overlay_path = "data/Zoning/dallas_historic_overlay/HistoricOverlay.shp"
dallas_historic_subdistricts_path = "data/Zoning/dallas_historic_subdistricts/HistoricSubdistricts.shp"

houston_historic_districts_path = "data/Zoning/houston_historic_protections/HISTORIC_DISTRICTS_CITY.shp"
houston_historic_landmarks_path = "data/Zoning/houston_historic_protections/HISTORICAL_SITES.shp"

# Houston special setbacks and lots
houston_spec_setbacks_path = "data/Zoning/houston_spec_minimum_building_lines/Special_Minimum_Building_Lines.shp"
houston_spec_min_lot_path = "data/Zoning/houston_spec_minimum_lot/Minimum_Lot_Size.shp"

# Civic clubs
houston_civic_clubs = "data/Zoning/houston_civic_clubs/cohgis.PD.shp"