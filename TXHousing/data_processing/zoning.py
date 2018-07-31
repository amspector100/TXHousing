"""Processes zoning data for Austin, Dallas, and all of North Texas"""

import time
import pandas as pd
import geopandas as gpd
import warnings

class zoning_inputs():

    """ Simple class which binds information about zoning datasets. Also, the lat/long features define the city
        center for each city."""

    def __init__(self, path, feature, separator, proj4string, base_zones, lat = 0, long = 0, zoom = 10,
                 title = '', crs = None, regulations_path = None):
        self.path = path
        self.feature = feature
        self.separator = separator
        self.proj4string = proj4string
        self.base_zones = base_zones
        self.lat = lat
        self.long = long
        self.zoom = zoom
        self.title = title
        self.crs = crs
        self.regulations_path = regulations_path

    def process_zoning_shapefile(self, overlay = {'-MU':'Multifamily'}, broaden = True, parse_base_zones = True, quietly = False):
        """
        :param broaden: Default true. If true, will decode zoning data into a broader classification (i.e. sf, mf)
            as specified by the zoning input class - this processed zoning data will go in a column labelled "broad_zone."
        :type broaden: Boolean
        :param overlay: A dictionary which maps overlay strings to broadened outputs. This is mostly useful for dealing with
            mixed use overlays in Austin (there are no mixed-use base zones, so passing in '-mu':'Multifamily'} will help
            the function recognize that mixed-use zones are multifamily zones.
        :param parse_base_zones: Boolean. Default true. If true, will use the regulations path provided by the input to
            parse the base zones of the zoning data.
        :param quietly: Default false. If true, will suppress all print statements and make the function "quiet"
        :type quietly: Boolean
        """

        # Begin timing
        time0 = time.time()

        # Read data and process zone codes
        if quietly is not True:
            print('Reading file')
        raw_data = gpd.read_file(self.path)
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

        if parse_base_zones:

            if self.regulations_path is not None:

                # Get regulations data - we only read this in because we need the list of base zones.
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

                # Apply
                raw_data['base_zone'] = raw_data[self.feature].apply(process_zone)

            else:
                warnings.warn(
                    """parse_base_zones = True, but a regulations path is required to parse base_zones from zoning data.""")

        # Manually set CRS if necessary, then transform
        if raw_data.crs is not None:
            raw_data.crs = self.crs
        raw_data = raw_data.to_crs({'init': 'epsg:4326'})

        if quietly is not True:
            print('Finished processing zones, mixed use overlay, and crs, took {}'.format(time.time() - time0))

        return raw_data

# Downloaded from http://data-nctcoggis.opendata.arcgis.com/datasets/2015-land-use
north_texas_inputs = zoning_inputs(path = 'data/Zoning Shapefiles/2015_North_Texas_Land_Use/2015_Land_Use.shp',
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
                                   title = 'Base Zones in and around Dallas, Texas')

# Downloaded from https://data.austintexas.gov/Locations-and-Maps/Zoning/5rzy-nm5e
austin_regulations_path = "shared_data/austin zoning standards.csv"
austin_inputs = zoning_inputs(path = "data/Zoning Shapefiles/austin_zoning/geo_export_571668ee-52f1-4ac9-a4e0-3b8bb348eae7.shp",
                              feature = 'zoning_zty',
                              separator = '-',
                              proj4string = 'EPSG:4326',
                              base_zones = {'Single Family':['SF'], 'Multifamily':['MF'],
                                            'Other Residential':['MH', 'RR', 'LA'], 'Other':[]},
                              lat = 30.267,
                              long = -97.743,
                              zoom = 9,
                              title = 'Base Zones in Austin, Texas',
                              crs = {'init': 'epsg:4326'},
                              regulations_path=austin_regulations_path)

# From dallas's open data site
dallas_regulations_path = 'shared_data/dallas zoning standards.csv'
dallas_inputs = zoning_inputs(path = "data/Zoning Shapefiles/DallasBaseZoning/BaseZoning.shp",
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
                              title = 'Base Zones in Dallas, Texas',
                              crs = {'init':'epsg:2276'},
                              regulations_path=dallas_regulations_path)

# Houston has no base zones but this is useful for making maps
houston_inputs = zoning_inputs(path = None,
                               feature = None,
                               separator = None,
                               proj4string = None,
                               base_zones = None,
                               lat = 29.7604,
                               long = -95.3698,
                               zoom = 9,
                               title = 'Houston Texas')