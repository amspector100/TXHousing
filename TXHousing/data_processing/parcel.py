import time
import pandas as pd
import geopandas as gpd
import re
import warnings
import os
from ..utilities import simple, spatial_joins, measurements
from . import boundaries, zoning

# Globals for municipalities -----------------------------------------------------------------------------------------
# Austin
austin_parcel_path = "data/Parcel/austin_parcels_2016/geo_export_813e97e4-7fde-4e3a-81b3-7ca9e8a89bd0.shp"

# Dallas
dallas_county_parcel_path = "data/Parcel/dalllas_county_2018/PARCEL/PARCEL.shp"# 2018
dallas_county_appraisal_path = 'data/Parcel/dallas_county_parcel_data/ACCOUNT_APPRL_YEAR.csv'
dallas_county_res_path = 'data/Parcel/dallas_county_parcel_data/res_detail.csv'
dallas_county_land_path = "data/Parcel/dallas_county_parcel_data/land.csv"

# Harris
harris_parcel_path_2018 = "data/Parcel/Harris_Parcels_2018/Parcels.shp"

# Extra parcel feature files and their columns for Harris - unfortunately the columns are not part of the data.
harris_parcel_land_path_2018 = "data/Parcel/Harris_Parcel_Land_Features/land.txt"
harris_parcel_appraisal_path_2018 = "data/Parcel/Harris_Parcel_Land_Features/real_acct.txt"
harris_parcel_building_res_path_2018 = 'data/Parcel/Harris_Parcel_Land_Features/building_res.txt'
harris_parcel_building_other_path_2018 = 'data/Parcel/Harris_Parcel_Land_Features/building_other.txt'
houston_land_columns = ["ACCOUNT", "LINE_NUMBER", "LAND_USE_CODE", "LAND_USE_DSCR", "SITE_CD", "SITE_CD_DSCR",
                        "SITE_ADJ", "UNIT_TYPE", "UNITS", "SIZE_FACTOR", "SITE_FACT", "APPR_OVERRIDE_FACTOR",
                        "APPR_OVERRIDE_REASON", "TOT_ADJ", "UNIT_PRICE", "ADJ_UNIT_PRICE", "VALUE", "OVERRIDE_VALUE"]
houston_appraisal_columns = ['ACCOUNT', 'TAX_YEAR', 'MAILTO', 'MAIL_ADDR_1', 'MAIL_ADDR_2',
                            'MAIL_CITY', 'MAIL_STATE', 'MAIL_ZIP', 'MAIL_COUNTRY', 'UNDELIV', 'STR_PFX',
                            'STR_NUM', 'STR_NUM_SFX', 'STR_NAME', 'STR_SFX', 'STR_SFX_DIR', 'STR_UNIT',
                            'SITE_ADDR_1', 'SITE_ADDR_2', 'SITE_ADDR_3', 'STATE_CLASS', 'SCHOOL_DIST',
                            'MAP_FACET', 'KEY_MAP', 'NEIGHBORHOOD_CODE', 'NEIGHBORHOOD_GROUP', 'MARKET_AREA_1',
                            'MARKET_AREA_1_DSCR', 'MARKET_AREA_2', 'MARKET_AREA_2_DSCR', 'ECON_AREA',
                            'ECON_BLD_CLASS', 'CENTER_CODE', 'YR_IMPR', 'YR_ANNEXED', 'SPLT_DT', 'DSC_CD',
                            'NXT_BUILDING', 'TOTAL_BUILDING_AREA', 'TOTAL_LAND_AREA', 'ACREAGE', 'CAP_ACCOUNT',
                            'SHARED_CAD_CODE', 'LAND_VALUE', 'IMPROVEMENT_VALUE', 'EXTRA_FEATURES_VALUE',
                            'AG_VALUE', 'ASSESSED_VALUE', 'TOTAL_APPRAISED_VALUE', 'TOTAL_MARKET_VALUE', 'PRIOR_LND_VALUE',
                             'PRIOR_IMPR_VALUE', 'PRIOR_X_FEATURES_VALUE', 'PRIOR_AG_VALUE', 'PRIOR_TOTAL_APPRAISED_VALLUE',
                            'PRIOR_TOTAL_MARKET_VALUE', 'NEW_CONSTRUCTION_VALUE', 'TOTAL_RCN_VALUE', 'VALUE_STATUS',
                            'NOTICED', 'NOTICE_DATE', 'PROTESTD', 'CERTIFIED_DATE', 'LAST_INSPECTED_DATE',
                            'LAST_INSPECTED_BY', 'NEW_OWNER_DATE', 'LEGAL_DSCR_1', 'LEGAL_DSCR_2', 'LEGAL_DSCR_3',
                            'LEGAL_DSCR_4', 'JURS']
houston_building_res_columns = ['ACCOUNT', 'USE_CODE', 'BUILDING_NUMBER', 'IMPRV_TYPE', 'BUILDING_STYLE_CODE',
                                'CLASS_STRUCTURE', 'CLASS_STRUC_DESCRIPTION', 'DEPRECIATION_VALUE',
                                'CAMA_REPLACEMENT_COST', 'ACCRUED_DEPR_PCT', 'QUALITY', 'QUALITY_DESCRIPTION',
                                'DATE_ERECTED', 'EFFECTIVE_DATE', 'YR_REMODEL', 'YR_ROLL', 'APPRAISED_BY',
                                'APPRAISED_DATE', 'NOTE', 'IMPR_SQ_FT', 'ACTUAL_AREA', 'HEAT_AREA', 'GROSS_AREA',
                                'EFFECTIVE_AREA', 'BASE_AREA', 'PERIMETER', 'PERCENT_COMPLETE', 'NBHD_FACTOR',
                                'RCNLD', 'SIZE_INDEX', 'LUMP_SUM_ADJ']
houston_building_other_columns = ["ACCOUNT", "USE_CODE", "BLD_NUM", "IMPRV_TYPE", "BUILDING_STYLE_CODE", "CLASS_STRUCTURE",
                                "CLASS_STRUC_DESCRIPTION", "NOTICED_DEPR_VALUE", "DEPRECIATION_VALUE", "MS_REPLACEMENT_COST",
                                "CAMA_REPLACEMENT_COST", "ACCRUED_DEPR_PCT", "QUALITY", "QUALITY_DESCRIPTION",
                                "DATE_ERECTED", "EFFECTIVE_DATE", "YR_REMODEL", "YR_ROLL", "APPRAISED_BY",
                                "APPRAISED_DATE", "NOTE", "IMPR_SQ_FT", "ACTUAL_AREA", "HEAT_AREA", "GROSS_AREA",
                                "EFFECTIVE_AREA", "BASE_AREA", "PERIMETER", "PERCENT_COMPLETE", "CATEGORY",
                                "CATEGORY_DSCR", "PROPERTY_NAME", "UNITS", "NET_RENT_AREA", "LEASE_RATE",
                                "OCCUPANCY_RATE", "TOTAL_INCOME"]

class Parcel(boundaries.Boundaries):
    """ Class for Parcel Data.

    :param path: The path of the data.
    :param account_col: The name of a column containing unique ids for each property.
    :param county: The county that the parcel data is in. Can also be a column of the data which lists the county of
        each parcel. Defaults to None but this param is highly recommended.
    :param processing_function: A custom processing function used to set as self.data
    :param geometry_column: The geometry column of the data, defaults to 'geometry'
    :param crs: If necessary, set the crs
    :param name: A name for the entire dataset, used to identify when printing. Defaults to None, at which point it
        uses the path of the dataset as its name (if this is also None then the name is just None).
    :param kwargs: kwargs to pass to the processing_function

    """

    def __init__(self, path, account_col, county = None, processing_function = None, geometry_column = 'geometry',
                 crs = None, name = None, **kwargs):

        # Record geometry column and initialization time
        self.time0 = time.time()
        self.geometry_column = geometry_column
        if name is not None:
            self.name = name
        else:
            self.name = path
        print('Staring to initialize at time {} for parcel {}'.format(time.time() - self.time0, self.name))


        # Custom processing function
        if processing_function is not None:
            self.data = processing_function(**kwargs)
        else:
            self.data = gpd.read_file(path)

        # Deal with crs
        if self.data.crs is None and crs is not None:
            self.data.crs = crs
        elif self.data.crs is None:
            warnings.warn('No crs set for this instance of parcel data')

        # Ensure validity of geometry
        self.data = self.data.loc[self.data[geometry_column].apply(lambda x: x is not None)]
        self.data = self.data.loc[self.data[geometry_column].is_valid]

        # Get the county of the data
        # County of the data
        if county is not None and county in self.data.columns:
            self.data = self.data.rename(columns = {county: 'county'})
        elif county is not None:
            self.data['county'] = county

        # Prop id column
        try:
            self.data['account'] = self.data[account_col].apply(lambda x: str(int(float((x)))))
        except ValueError:
            self.data = self.data.rename(columns = {account_col: 'account'})

        # Drop duplicates by centroids and also by account number. Centroid operations are very cheap so it's actually
        # more efficient to do this before transforming to lat long and then dropping centroids/centroids string and
        # doing it again.
        self.data = self.data.drop_duplicates(subset = 'account', keep = 'first')

        if 'centroids' not in self.data.columns:
            self.data['centroids'] = self.data[geometry_column].centroid
        self.data['centroids_string'] = self.data['centroids'].astype(str)
        self.data = self.data.drop_duplicates(subset = ['centroids_string'], keep = 'first')
        self.data = self.data.drop(['centroids', 'centroids_string'], axis = 1)
        print('Finished initializing at time {} for parcel {}'.format(time.time() - self.time0, self.name))

    def parse_broad_zone(self, description_column, broad_zone_dictionary):
        """
        Parses broad zones of parcel data. Updates the self.data attribute to include a 'broad_zone' column.

        :param description_column: The column of descriptions to parse
        :param broad_zone_dictionary: A dictionary mapping base zones (i.e. 'Single Family') to a list of strings
            which signal that a description means that base zone. Ex: {'Single Family':['sf', 'single f'], 'Multifamily':['mf']}
            Note that the order of this dictionary DOES matter - the function will return the FIRST key which has a match.
        """

        self.data['zone_feature'] = self.data[description_column].fillna('').astype(str)

        # Create helper function
        def process_broad_zone(text):
            for key in broad_zone_dictionary:
                # This checks whether any of the strings in the broad_zone_dictionary[key] list appear in the text
                if bool(re.search(('|').join(broad_zone_dictionary[key]), text)):
                    return key
            return "Other"

        # Apply
        self.data['broad_zone'] = self.data['zone_feature'].apply(process_broad_zone)

    def merge_multiple(self, merge_paths, right_keys, left_keys = None, **kwargs):
        """A wrapper of the geopandas.merge function, quickly merge with multiple other csvs

        :param merge_paths: If the geodata must be merged with another datasource, this should be an iterable containing
            the aths of the data to merge it with.
        :param left_keys: If the geodata must be merged with another datasource, this should be an iterable containing the
            left keys used to merge the data, in the same order of as the 'merge_paths'. Defaults to None, in which case
            the merge is performed using the account_feature as the merge key.
        :param right_keys: If the geodata must be merged with another datasource, this should be an iterable containing the
            right keys used to merge the data, in the same order of as the 'merge_paths'
        :param kwargs: kwargs to pass to the .merge call

        """

        #  Start by ensuring that the left_keys, right_keys, merge_paths are iterables
        if isinstance(merge_paths, str):
            merge_paths = [merge_paths]
        if isinstance(left_keys, str):
            left_keys = [left_keys]
        elif left_keys is None:
            left_keys = ['account'] * len(merge_paths)
        if isinstance(right_keys, str):
            right_keys = [right_keys]
        if len(merge_paths) != len(right_keys) or len(merge_paths) != len(left_keys):
            raise ValueError('The lengths of left_keys ({}), right_keys ({}), and merge_paths ({}) must be equal'.format(
                len(left_keys),
                len(right_keys),
                len(merge_paths)))

        # Read data and merge
        for merge_path, left_key, right_key in zip(merge_paths, left_keys, right_keys):
            extra_data = pd.read_csv(merge_path, **kwargs)
            extra_data = extra_data.loc[(extra_data[right_key].notnull())]
            extra_data = extra_data.drop_duplicates(subset=right_key, keep='first')
            self.data = self.data.merge(extra_data, left_on=left_key, right_on=right_key, how='left')

    def measure_parcels(self, lat, long, area_feature = None):
        """ Calculates the distance to center of the city as well as the area and centroid of each parcel. In effect,
        it adds four new columns to self.data: lat, long, centroids, area_sqft, and dist_to_center. This also
        transforms the parcels to lat long.

        """

        # If we can, just use the area feature
        if area_feature is not None:
            self.data = self.data.rename(columns = {area_feature: 'area_sqft'})

            if self.data.crs != {'init':'epsg:4326'}:
                self.data = self.data.to_crs({'init':'epsg:4326'})

        # Else we have to transform to lat long.
        else:
            self.data = measurements.get_area_in_units(self.data, scale = 1, name = 'area_sqft',
                                                       final_projection = {'init':'epsg:4326'})

        # This generates centroids, lat, long, and dist to center.
        self.data['dist_to_center'] = measurements.calculate_dist_to_center(self.data, lat, long, drop_centroids=False)

    def pull_geographic_information(self, bounding_counties):

        # Actually calculate spatial index for efficiency
        self.data = simple.process_geometry(self.data, geometry_column = self.geometry_column)
        if 'centroids' not in self.data.index:
            self.data['centroids'] = self.data[self.geometry_column].centroid
        self.data = self.data.set_geometry('centroids')
        self.centroid_spatial_index = self.data.sindex
        self.data = self.data.set_geometry(self.geometry_column)

        # Pull ua, and place information
        for path, colname, fill in zip([boundaries.ua_path, boundaries.texas_places_path],
                                       ['ua', 'place'],
                                       ['None', 'Unincorporated']):

            # Automatically converts to lat/long
            boundary_class = boundaries.Boundaries(path, bounding_counties = bounding_counties, to_latlong = True)
            boundary_class.data = boundary_class.data.rename(columns = {'NAME10':'NAME'})

            # Get intersections
            print('Taking intersections for {}'.format(colname))
            mapper = boundary_class.fast_intersection(self.data, geometry_column = self.geometry_column,
                                                      small_points_spatial_index = self.centroid_spatial_index,
                                                      horiz = 20, vert = 20)
            self.data[colname] = self.data.index.map(mapper)
            self.data[colname] = self.data[colname].map(boundary_class.data['NAME']).fillna(fill).astype(str)


        # Do zip separately because it has its own class
        zipdata = boundaries.ZipBoundaries(bounding_counties = bounding_counties, to_latlong = True)
        print('Taking intersections for zipcode')
        mapper = zipdata.fast_intersection(self.data, geometry_column=self.geometry_column,
                                           small_points_spatial_index=self.centroid_spatial_index,
                                           horiz=20, vert=20)
        self.data['zipcode'] = self.data.index.map(mapper)
        self.data['zipcode'] = self.data['zipcode'].fillna('Other').astype(str)

    def process_parcel_data(self, broad_zone_feature, broad_zone_dictionary, zoning_input, bounding_counties,
                            area_feature = None, merge_paths = None, left_keys = None, right_keys = None,
                            save_path = None, geo_save_path = None):
        """Wrapper which calls self.merge_multiple, self.parse_broad_zone, self.measure_parcels,
        and self.pull_geographic_information in that order. Basically, after initializing the parcel data and calling
        this function, it should have the following features:

        account, broad_zone, zone_feature, dist_to_center, area_sqft, lat, long, county, zipcode, place, ua (urban area)

        as well as a host of other features that may have been joined with/initially part of the data.

        :param broad_zone_feature: The feature used to parse broad_zones. State cds codes are preferred for consistency.
        :param broad_zone_dictionary: Dictionary mapping broad zones to keywords that we will use to parse the
            broad_zone_feature and obtain broad zones.
        :param zoning_input: The zoning input for the city of interest; used to calculate distance from city center.
        :param bounding_counties: A list of counties which might intersect the data. It's computationally relatively
            inexpensive to add extra counties to this list, by the way.
        :param area_feature: Default None. If the data already lists its area in square feet, then we won't bother
            recalculating it (this saves a ton of time because crs transformations are very expensive for parcel data).
        :param merge_paths: If the geodata must be merged with another datasource, this should be an iterable containing the
            paths of the data to merge it with.
        :param left_keys: If the geodata must be merged with another datasource, this should be an iterable containing the
            left keys used to merge the data, in the same order of as the 'merge_paths'. Defaults to None, in which case
             the merge is performed using 'account', which is parsed by the init function, as the merge key.
        :param right_keys: If the geodata must be merged with another datasource, this should be an iterable containing the
            right keys used to merge the data, in the same order of as the 'merge_paths'
        :param save_path: A csv path at which to save the data. Defaults to None, in which case it will not save the
            data.
        :param geo_save_path: A .shp path at which to save the data. Defaults to None, in which case it will not save
            the data.
        :return: None, but modifies parcel data in place.
        """

        # Merge with external data
        print('Starting to merge and parse at time {} for {} parcels'.format(time.time() - self.time0, self.name))
        if merge_paths is not None:
            self.merge_multiple(merge_paths = merge_paths, left_keys = left_keys, right_keys = right_keys)

        # Process broad_zone
        self.parse_broad_zone(description_column = broad_zone_feature, broad_zone_dictionary = broad_zone_dictionary)

        # Get centroids, dist to center, lat, long, etc
        print('Starting to measure distances and areas at time {} for {} parcels'.format(time.time() - self.time0, self.name))
        self.measure_parcels(lat = zoning_input.lat, long = zoning_input.long, area_feature = area_feature)

        # Get geographic information
        print('Starting to pull geographic info at time {} for {} parcels'.format(time.time() - self.time0, self.name))
        self.pull_geographic_information(bounding_counties)

        # Save to csv if save_path is not None:
        if save_path is not None:
            print('Saving parcel data to csv at time {} for {} parcels'.format(time.time() - self.time0, self.name))
            csv_data = self.data[[c for c in self.data.columns if c not in ['centroids', 'geometry']]]
            csv_data.to_csv(save_path)

        # Save to .shp if geo_save_path is not None:
        if geo_save_path is not None:
            print('Saving parcel data to .shp at time {} for {} parcels'.format(time.time() - self.time0, self.name))
            shp_data = self.data[[c for c in self.data.columns if c not in ['centroids']]]
            shp_data.to_file(geo_save_path)

        # Print success!
        print('Finished processing {} parcel shapefile; took {} since initialization'.format(self.name,
                                                                                             time.time() - self.time0))

# Specific processing functions --------------------------------------------------------------------------------------

def process_austin_parcel_data():
    """Reads Austin parcel data and processes base zones."""

    time0 = time.time()
    print('Starting to read parcel data')
    parcel_data = gpd.read_file(austin_parcel_path)
    print('Finished reading parcel data, took {} seconds. Now processing base zones.'.format(time.time() - time0))

    # Get basezone in proper format - originally it provides a description in parenthesis as well which is not necessary
    def format_basezone(s):
        s = str(s)
        try:
            return s.split('(')[1][:-1]
        except:
            # Usually this means there's no basezone information
            return s
    parcel_data['basezone'] = parcel_data['basezone'].apply(format_basezone)

    return parcel_data

def process_houston_parcel_data(feature_files = [harris_parcel_building_res_path_2018],
                                feature_columns_list = [houston_building_res_columns],
                                process_crs = True,
                                county_level = False):
    """ Merge houston or harris parcel data with harris county data.

    :param feature_files: A list of paths of feature files to merge with the data. Defaults to
        the building_res file path. Can also be a string of one file (doesn't have to be a list). If you set this equal
        to None, it will just return the parcel shape data with no other data attached.
    :param feature_columns_list: A list of column headers for each file - you can find these in the
        Harris_Parcel_Feature_Columns microsoft database under data.
    :param bool process_crs: Default True. If True, transform parcel data to latlong. This is very expensive.
    :param county_level: Default False. If False, subset the data to only consider the parcels within the Houston
        municipal boundaries.
    :return: GeoDataFrame with the parcel shapes and the merged data.

    """

    # Get parcel data
    geodata = gpd.read_file(harris_parcel_path_2018)
    geodata = geodata.loc[geodata['geometry'].is_valid]
    geodata.reset_index(inplace = True)
    geodata = geodata.loc[geodata['HCAD_NUM'].apply(simple.will_it_float)]
    geodata['HCAD_NUM'] = geodata['HCAD_NUM'].astype(float)

    # Make sure feature_files is a list not a string
    if isinstance(feature_files, str):
        feature_files = [feature_files]

    # Get data
    if feature_files is not None:
        for feature_path, feature_header in zip(feature_files, feature_columns_list):
            data = pd.read_csv(feature_path, sep = '\t', header = None, encoding='Latin-1')
            data.columns = feature_header
            data = data.rename(columns = {'ACCOUNT':'HCAD_NUM'})
            data = data.drop_duplicates(subset=['HCAD_NUM'], keep='first')
            data['HCAD_NUM'] = data['HCAD_NUM'].astype(float)
            geodata = geodata.merge(data, how = "left", on = "HCAD_NUM")

    # Subset to only include houston if indicated
    if county_level == False:
        place_shapes = gpd.read_file(boundaries.texas_places_path)
        place_shapes = gpd.GeoDataFrame(geometry = place_shapes.loc[place_shapes['NAME'] == 'Houston'].to_crs(geodata.crs).simplify(tolerance = 0.005))
        mapper = spatial_joins.fast_polygon_intersection(geodata, place_shapes, horiz = 25, vert = 25)
        geodata['place'] = geodata.index.map(mapper)
        geodata = geodata.loc[geodata['place'].notnull()]

    # Transform to lat long
    if process_crs:
        geodata = geodata.to_crs({'init':'epsg:4326'})

    # Return
    return geodata

# Cache parcel data --------------------------------------------------------------------------------------------------

def get_cached_municipal_parcel_path(cityname):
    return 'data/Caches/municipal_parcel/{}_municipal_parcel.csv'.format(cityname)

def get_cached_all_parcel_path_csv(cityname):
    return 'data/Caches/all_parcel/{}_all_parcel.csv'.format(cityname)

def get_cached_all_parcel_path_shp(cityname):
    return 'data/Caches/all_parcel/{}_all_parcel.shp'.format(cityname)

all_dallas_zoning_path_csv = 'data/Zoning/processed_north_texas_data/zones.csv'

# State code dictionary is used for a lot of different counties, mostly around Houston
state_sptbcode_dictionary = {'Single Family':['A1', 'A2'], 'Multifamily':['A3', 'A4', 'B1', 'B2', 'B3', 'B4']}

# Create caches of parcel data to make reading/working with it easier
def cache_municipal_parcel_data():
    """Creates csvs which store the centroids, area, and other relevant features about each parcel for municipalities."""
    time0 = time.time()

    # Calling the processing functions, then saving as csvs.


    # Houston -----
    houston_parcels = Parcel(path = None,
                             account_col = 'HCAD_NUM',
                             county = 'Harris',
                             name = 'Houston',
                             processing_function = process_houston_parcel_data,
                             feature_files=[harris_parcel_land_path_2018,
                                            harris_parcel_appraisal_path_2018,
                                            harris_parcel_building_res_path_2018],
                             feature_columns_list=[houston_land_columns, houston_appraisal_columns,
                                                   houston_building_res_columns],
                             county_level=False)

    houston_parcels.process_parcel_data(broad_zone_feature = 'STATE_CLASS',
                                       broad_zone_dictionary = state_sptbcode_dictionary,
                                       zoning_input = zoning.houston_inputs,
                                       bounding_counties = ['Harris'],
                                       area_feature = None,
                                       merge_paths = None,
                                       left_keys = None,
                                       right_keys = None,
                                       save_path = get_cached_municipal_parcel_path('houston'))

    # Austin -----
    print('Starting to process Austin parcel data at global time {}'.format(time.time() - time0))
    # Read data and process base zone
    austin_parcels = Parcel(path = None,
                            account_col = 'prop_id',
                            county = 'Travis',
                            name = 'Austin',
                            processing_function = process_austin_parcel_data)

    # Note SF-4B refers to condominium, and SF-6 refers to 'townhouse and condominium'
    austin_dictionary = {'Single Family': ['SF-1', 'SF-2', 'SF-3', 'SF-4A', 'SF-5'],
                         'Multifamily': ['SF-4B', 'SF-6', 'MF-1', 'MF-2', 'MF-3', 'MF-4', 'MF-5', 'MF-6']}
    austin_parcels.process_parcel_data(broad_zone_feature = 'basezone',
                                       broad_zone_dictionary = austin_dictionary,
                                       zoning_input = zoning.austin_inputs,
                                       bounding_counties = ['Travis'],
                                       area_feature='shape_area',
                                       merge_paths = None,
                                       left_keys = None,
                                       right_keys = None,
                                       save_path = get_cached_municipal_parcel_path('austin'))

    # Dallas --
    print('Starting to process Dallas parcel data at global time {}'.format(time.time() - time0))
    dallas_parcels = Parcel(path = dallas_county_parcel_path, account_col = 'Acct', county = 'Dallas', name = 'Dallas')
    dallas_sptb_dictionary = {'Single Family': ['A11'], 'Multifamily': ['B11', 'B12', 'A12', 'A13']}
    dallas_parcels.process_parcel_data(broad_zone_feature = 'SPTD_CODE',
                                       broad_zone_dictionary = dallas_sptb_dictionary,
                                       zoning_input = zoning.dallas_inputs,
                                       bounding_counties = ['Dallas'],
                                       area_feature='Shape_area',
                                       merge_paths=[dallas_county_appraisal_path, dallas_county_res_path],
                                       left_keys = None,
                                       right_keys=['ACCOUNT_NUM', 'ACCOUNT_NUM'],
                                       save_path = get_cached_municipal_parcel_path('dallas'))

def cache_all_parcel_data():
    """Creates csvs which store the centroids, area, and other relevant features about each parcel for all of the
    surrounding counties of each core municipality."""


    time0 = time.time()
    final_columns = ['account', 'broad_zone', 'zone_feature', 'dist_to_center', 'area_sqft', 'lat', 'long', 'county',
                     'zipcode', 'place', 'ua', 'geometry']

    # Houston  --------------------------------------------------------------------------------------------
    # All of these counties use the state code dictionary - manually checked that this is the right thing to do

    fort_bend_parcel_path = 'data/Parcel/fort_bend_parcels_2018/CAMASUMMARY.shp'
    if os.path.exists(fort_bend_parcel_path) == False:
        raise FileNotFoundError("""Fort Bend County Parcel not in the data directory - use cached data instead or 
            download the raw data from  'https://fbcad.org/District-Information/GIS-Data'.""")

    montgomery_county_parcel_path = "data/Parcel/montgomery_parcels_2018/Tax_Parcel_View.shp"  # 2018
    if os.path.exists(montgomery_county_parcel_path) == False:
        raise FileNotFoundError("""Montgomery County Parcel not in the data directory - use cached data instead or 
            download the raw data from 'https://data-moco.opendata.arcgis.com/datasets/tax-parcel-view'.""")

    fort_bend_parcels = Parcel(path = fort_bend_parcel_path, account_col = 'UID',
                               county = 'Fort Bend', name = 'Fort Bend')
    fort_bend_parcels.process_parcel_data(broad_zone_feature = 'LMainSegSP',
                                          broad_zone_dictionary = state_sptbcode_dictionary,
                                          zoning_input = zoning.houston_inputs,
                                          bounding_counties = ['Fort Bend'],
                                          area_feature = 'Shape_Area',
                                          merge_paths = None,
                                          left_keys = None,
                                          right_keys = None)
    fort_bend_parcels.data = fort_bend_parcels.data[final_columns]

    montgomery_parcels = Parcel(path = montgomery_county_parcel_path, account_col = 'PropertyNu',
                                county = "Montgomery", name = 'Montgomery')
    montgomery_parcels.process_parcel_data(broad_zone_feature = 'fStateCode',
                                           broad_zone_dictionary = state_sptbcode_dictionary,
                                           zoning_input = zoning.houston_inputs,
                                           bounding_counties = ['Montgomery'],
                                           merge_paths = None)
    montgomery_parcels.data = montgomery_parcels.data[final_columns]



    # And lastly houston, which takes forever
    houston_parcels = Parcel(path = None,
                             account_col = 'HCAD_NUM',
                             county = 'Harris',
                             name = 'Houston',
                             processing_function = process_houston_parcel_data,
                             feature_files=[harris_parcel_land_path_2018,
                                            harris_parcel_appraisal_path_2018,
                                            harris_parcel_building_res_path_2018],
                             feature_columns_list=[houston_land_columns, houston_appraisal_columns,
                                                   houston_building_res_columns],
                             county_level=True)
    houston_parcels.process_parcel_data(broad_zone_feature = 'STATE_CLASS',
                                       broad_zone_dictionary = state_sptbcode_dictionary,
                                       zoning_input = zoning.houston_inputs,
                                       bounding_counties = ['Harris'],
                                       area_feature = None,
                                       merge_paths = None,
                                       left_keys = None,
                                       right_keys = None)
    houston_parcels.data = houston_parcels.data[final_columns]

    print('Combining and saving for Houston, time is {}'.format(time.time() - time0))
    all_houston_parcels = pd.concat([houston_parcels.data, montgomery_parcels.data, fort_bend_parcels.data],
                                   axis=0,
                                   ignore_index=True)
    all_houston_parcels.to_file(get_cached_all_parcel_path_shp('houston'))
    all_houston_parcels_csv = all_houston_parcels[[col for col in all_houston_parcels.columns if col != 'geometry']]
    all_houston_parcels_csv.to_csv(get_cached_all_parcel_path_csv('houston'))
    print('Finished with Houston, global time is {}'.format(time.time() - time0))

    # Austin  --------------------------------------------------------------------------------------------
    travis_county_parcel_path = "data/Parcel/travis_county_parcels_2016/Parcels_Travis_2016.shp"  # 2016
    travis_county_data_path = 'data/Parcel/travis_county_parcel_data/land_det.csv'
    if (os.path.exists(travis_county_parcel_path) and os.path.exists(travis_county_data_path)) == False:
        raise FileNotFoundError("""Travis County Parcel Shapes/Data not in the data directory - use cached data 
            instead or download the raw data from https://www.traviscad.org/reports-request/""")

    williamson_county_parcel_path = "data/Parcel/williamson_parcels_2016/Parcel_Poly.shp"  # 2017
    williamson_county_real_improvement_path = "data/Parcel/williamson_data_2016b/Improvement.txt"  # 2018
    if (os.path.exists(williamson_county_parcel_path) and
        os.path.exists(williamson_county_real_improvement_path)) == False:
        raise FileNotFoundError("""Williamson County Parcel Shapes/Data not in the data directory - use cached data 
            instead or download the raw data from https://www.wcad.org/data-downloads/""")


    williamson_parcels = Parcel(path = williamson_county_parcel_path, account_col = 'PIN',
                                county = 'Williamson', name = 'Williamson')
    # See https://www.wcad.org/wp-content/uploads/2016/09/2015Report.pdf for the dictionary information
    williamson_dictionary = {'Single Family':['A1', 'A9'], 'Multifamily':['A8', 'B1', 'B2', 'B4']}
    williamson_parcels.process_parcel_data(broad_zone_feature = 'StateCode',
                                           broad_zone_dictionary = williamson_dictionary,
                                           zoning_input = zoning.austin_inputs,
                                           bounding_counties = ['Williamson'],
                                           area_feature='Shape_area',
                                           merge_paths=[williamson_county_real_improvement_path],
                                           left_keys=None,
                                           right_keys=['QuickRefID'])
    williamson_parcels.data = williamson_parcels.data[final_columns]

    # See https://tax-office.traviscountytx.gov/pages/SPTC.php - A4/A5 refer to condos.
    travis_parcels = Parcel(path = travis_county_parcel_path, account_col = 'PROP_ID',
                            county = 'Travis', name = 'Travis')
    travis_dictionary = {'Single Family':['A1', 'A2', 'A3'], 'Multifamily':['A4', 'A5', 'B1', 'B2', 'B3', 'B4']}
    travis_parcels.process_parcel_data(broad_zone_feature = 'state_cd',
                                       broad_zone_dictionary = travis_dictionary,
                                       zoning_input = zoning.austin_inputs,
                                       bounding_counties = ['Travis'],
                                       area_feature = 'Shape__Are',
                                       merge_paths=[travis_county_data_path],
                                       left_keys=['PROP_ID'],
                                       right_keys=['PROP_ID'])
    travis_parcels.data = travis_parcels.data[final_columns]

    print('Combining and saving for Austin, global time is {}'.format(time.time() - time0))
    all_austin_parcels = pd.concat([travis_parcels.data, williamson_parcels.data],
                                   axis=0,
                                   ignore_index=True)
    all_austin_parcels.to_file(get_cached_all_parcel_path_shp('austin'))
    all_austin_parcels_csv = all_austin_parcels[[col for col in all_austin_parcels.columns if col != 'geometry']]
    all_austin_parcels_csv.to_csv(get_cached_all_parcel_path_csv('austin'))
    print('Finished with Austin, global time is {}'.format(time.time() - time0))



    # Dallas -----------------------------------------------------------------------------------------------

    # Extra paths for the others. Make sure they exist or raise errors otherwise.
    collin_county_parcel_path = "data/Parcel/collin_county_2018/parcels.shp"  # 2018
    if os.path.exists(collin_county_parcel_path) == False:
        raise FileNotFoundError("""Collin County Parcel not in the data directory - use cached data instead or 
            download the raw data from https://www.collincad.org/downloads.""")

    denton_county_parcel_path = "data/Parcel/denton_county_parcels_2018/County_Parcels.shp"  # 2018
    if os.path.exists(denton_county_parcel_path) == False:
        raise FileNotFoundError("""Denton County Parcel not in the data directory - use cached data instead or 
            download the raw data from https://www.dentoncad.com/forms-and-downloads.""")

    processed_tarrant_county_parcel_path = "data/Parcel/processed_tarrant_county_parcels_2018/TADData.shp"
    if os.path.exists(processed_tarrant_county_parcel_path) == False:
        raise FileNotFoundError("""Processed Tarrant County Parcel data not in the data directory. Use the cached
            data or download the raw data from https://www.tad.org/data-download/ and run the 
            tarrant-parcel-processing.R script (this last step is necessary to avoid a bug in GeoPandas).""")

    print('Starting to process Tarrant parcel data at global time {}'.format(time.time() - time0))
    tarrant_parcels = Parcel(path = processed_tarrant_county_parcel_path,
                             account_col = 'TAXPIN', county = 'Tarrant', name = 'Tarrant')
    tarrant_parcels.process_parcel_data(broad_zone_feature = 'Prprt_C',
                                        broad_zone_dictionary = state_sptbcode_dictionary,
                                        zoning_input = zoning.dallas_inputs,
                                        bounding_counties = ['Tarrant'],
                                        area_feature='Lnd_SqF',
                                        merge_paths=None,
                                        left_keys = None,
                                        right_keys=None)
    tarrant_parcels.data = tarrant_parcels.data[final_columns]

    # A3 refers to Condos, A4 refers to townhomes. Note I have not included property improvements (classes A6, A9, B6, B9)
    # - see the data dictionary in the data/parcels directory.
    print('Starting to process Collin parcel data at global time {}'.format(time.time() - time0))
    collin_parcels = Parcel(path = collin_county_parcel_path, account_col = 'PROP_ID',
                            county = 'Collin', name = 'Collin')
    collin_dictionary = {'Single Family': ['A1'], 'Multifamily': ['A3', 'A4', 'B1', 'B2', 'B3', 'B4']}
    collin_parcels.process_parcel_data(broad_zone_feature = 'state_cd',
                                       broad_zone_dictionary = collin_dictionary,
                                       zoning_input = zoning.dallas_inputs,
                                       bounding_counties = ['Collin'],
                                       area_feature=None,
                                       merge_paths=None,
                                       left_keys = None,
                                       right_keys=None)
    collin_parcels.data = collin_parcels.data[final_columns]


    print('Starting to process Denton parcel data at global time {}'.format(time.time() - time0))
    denton_parcels = Parcel(path = denton_county_parcel_path, account_col = 'PROP_ID',
                            county = 'Denton', name = 'Denton')
    denton_dictionary = {'Single Family': ['Single Family', 'Duplexes'],
                         'Multifamily': ['Apartment', 'Townhomes', ', Condos'],
                         'Other Residential': ['Residential']}
    denton_parcels.process_parcel_data(broad_zone_feature = 'CD_DESCRIP',
                                       broad_zone_dictionary = denton_dictionary,
                                       zoning_input = zoning.dallas_inputs,
                                       bounding_counties = ['Denton'],
                                       area_feature = 'Shape__Are',
                                       merge_paths = None,
                                       left_keys = None,
                                       right_keys= None)
    denton_parcels.data = denton_parcels.data[final_columns]


    print('Starting to process Dallas parcel data at global time {}'.format(time.time() - time0))
    dallas_parcels = Parcel(path = dallas_county_parcel_path, account_col = 'Acct', county = 'Dallas', name = 'Dallas')
    dallas_sptb_dictionary = {'Single Family': ['A11'], 'Multifamily': ['B11', 'B12', 'A12', 'A13']}
    dallas_parcels.process_parcel_data(broad_zone_feature = 'SPTD_CODE',
                                       broad_zone_dictionary = dallas_sptb_dictionary,
                                       zoning_input = zoning.dallas_inputs,
                                       bounding_counties = ['Dallas'],
                                       area_feature='Shape_area',
                                       merge_paths=[dallas_county_appraisal_path, dallas_county_res_path],
                                       left_keys = None,
                                       right_keys=['ACCOUNT_NUM', 'ACCOUNT_NUM'])
    dallas_parcels.data = dallas_parcels.data[final_columns]

    print('Combining and saving for Dallas, time is {}'.format(time.time() - time0))
    all_dallas_parcels = pd.concat([collin_parcels.data, dallas_parcels.data, denton_parcels.data, tarrant_parcels.data],
                                   axis=0,
                                   ignore_index=True)
    all_dallas_parcels.to_file(get_cached_all_parcel_path_shp('dallas'))
    all_dallas_parcels_csv = all_dallas_parcels[[col for col in all_dallas_parcels.columns if col != 'geometry']]
    all_dallas_parcels_csv.to_csv(get_cached_all_parcel_path_csv('dallas'))
    print('Finished with Dallas, global time is {}'.format(time.time() - time0))

def cache_north_texas_zoning_data():
    """  Because we use the north texas zoning data to validate parcel results, we apply the parcel processing functions
     to the north texas zoning data and then cache the results.(This is a bit hacky, because you wouldn't expect to see
     zoning data fit under the parcel class, but so be it.) """

    final_columns = ['account', 'broad_zone', 'zone_feature', 'dist_to_center', 'area_sqft', 'lat', 'long', 'county',
                     'zipcode', 'place', 'ua']

    def retrieve_north_texas_data():
        north_texas_data = zoning.north_texas_inputs.process_zoning_shapefile(broaden = False, to_latlong = False)
        north_texas_data['id'] = north_texas_data.index # The data has no ID so we make a unique one
        return north_texas_data

    north_texas_data = Parcel(path = None, account_col = 'id', county = 'COUNTY', name = 'NT_Zoning',
                              processing_function = retrieve_north_texas_data)

    north_texas_data.process_parcel_data(broad_zone_feature = 'CATEGORY',
                                         broad_zone_dictionary  = zoning.north_texas_inputs.base_zones,
                                         zoning_input = zoning.dallas_inputs,
                                         bounding_counties = ['Erath', 'Johnson', 'Navarro', 'Hood', 'Palo Pinto',
                                                                        'Parker', 'Tarrant', 'Kaufman', 'Somervell', 'Wise',
                                                                        'Collin', 'Hunt', 'Rockwall', 'Dallas', 'Denton',
                                                                        'Ellis'],
                                         merge_paths = None)
    north_texas_data.data = north_texas_data.data[final_columns]

    north_texas_data.data.to_csv(all_dallas_zoning_path_csv)
