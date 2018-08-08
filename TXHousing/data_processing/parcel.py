import time
import pandas as pd
import geopandas as gpd
import re
import warnings
from ..utilities import simple, spatial_joins, measurements
from . import boundaries, zoning

# Globals ------------------------------------------------------------------------------------------------------------
austin_parcel_path = "data/Zoning Shapefiles/Austin Land Database 2016/geo_export_813e97e4-7fde-4e3a-81b3-7ca9e8a89bd0.shp"

dallas_parcel_data_path_2013 = "data/Zoning Shapefiles/Dallas County Parcels 2013/geo_export_9b090abf-d5d9-4c74-a6be-4486e75ee147.shp"
dallas_parcel_data_path_2016 = "data/Zoning Shapefiles/Dallas County Parcels 2016/geo_export_bd65f212-41ce-4166-804a-5dc5ef85ee84.shp"

harris_parcel_path_2018 = "data/Zoning Shapefiles/Harris_Parcels_2018/Parcels.shp"

# Extra parcel feature files and their columns for Harris - unfortunately these are not part of the data.
harris_parcel_land_path_2018 = "data/Zoning Shapefiles/Harris_Parcel_Land_Features/land.txt"
harris_parcel_appraisal_path_2018 = "data/Zoning Shapefiles/Harris_Parcel_Land_Features/real_acct.txt"
harris_parcel_building_res_path_2018 = 'data/Zoning Shapefiles/Harris_Parcel_Land_Features/building_res.txt'
harris_parcel_building_other_path_2018 = 'data/Zoning Shapefiles/Harris_Parcel_Land_Features/building_other.txt'

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

        self.data[description_column] = self.data[description_column].fillna('').astype(str)

        # Create helper function
        def process_broad_zone(text):
            for key in broad_zone_dictionary:
                # This checks whether any of the strings in the broad_zone_dictionary[key] list appear in the text
                if bool(re.search(('|').join(broad_zone_dictionary[key]), text)):
                    return key
            return "Other"

        # Apply
        self.data['broad_zone'] = self.data[description_column].apply(process_broad_zone)

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
        it adds six new columns to self.data: lat, long, centroids, area_sqft, and dist_to_center. This also
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
                            save_path = None):
        """Wrapper which calls self.merge_multiple, self.parse_broad_zone, self.measure_parcels,
        and self.pull_geographic_information in that order. Basically, after initializing the parcel data and calling
        this function, it should have the following features:
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

        # Save if save_path is not None:
        if save_path is not None:
            print('Saving parcel data at time {} for {} parcels'.format(time.time() - self.time0, self.name))
            csv_data = self.data[[c for c in self.data.columns if c not in ['centroids', 'geometry']]]
            csv_data.to_csv(save_path)


        # Print success!
        print('Finished processing {} parcel shapefile; took {} since initialization'.format(self.name,
                                                                                             time.time() - self.time0))


    # Basic processing functions -----------------------------------------------------------------------------------------
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

def join_dallas_parcel_data():
    """Joins 2013 and 2016 Dallas parcel data. Hopefully this will never be necessary as it takes a while."""

    time0 = time.time()
    print('Reading 2013 file')
    p2013 = gpd.read_file(dallas_parcel_data_path_2013)
    p2013.drop_duplicates('gis_acct', inplace = True) # Only drops about 700/300K
    p2013.columns = [str(col) + '_2013' for col in p2013.columns]
    p2013 = p2013.rename(columns = {'gis_acct_2013': 'gis_acct'})

    print('Finished reading 2013 file, took {}. Now reading 2016 file.'.format(time.time() - time0))
    p2016 = gpd.read_file(dallas_parcel_data_path_2016)
    p2016.drop_duplicates('gis_acct', inplace = True) # Only drops about 750/400K
    p2016.columns = [str(col) + '_2016' for col in p2016.columns]
    p2016 = p2016.rename(columns = {'gis_acct_2016': 'gis_acct'})

    print('Finished reading 2016 file, took {}. Now merging files.'.format(time.time() - time0))
    parcel_data = p2016.merge(p2013, on = 'gis_acct')
    parcel_data = gpd.GeoDataFrame(data = parcel_data[[col for col in parcel_data.columns if 'geometry' not in col]],
                                   geometry = parcel_data['geometry_2016'])
    print(parcel_data.columns, parcel_data.index)
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

# State code dictionary is used for a lot of different counties, mostly around Houston
state_sptbcode_dictionary = {'Single Family':['A1', 'A2'], 'Multifamily':['A3', 'A4', 'B1', 'B2', 'B3', 'B4']}

# Create caches of parcel data to make reading/working with it easier
def cache_municipal_parcel_data():
    """Creates csvs which store the centroids, area, and other relevant features about each parcel."""
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

def get_cached_municipal_parcel_path(cityname):
    return 'data/caches/municipal_parcel/{}_municipal_parcel.csv'.format(cityname)



dallas_county_parcel_path = "data/parcels/dalllas_county_2018/PARCEL/PARCEL.shp"# 2018
dallas_county_land_path = "data/parcels/dallas_county_parcel_data/land.csv"
dallas_county_appraisal_path = 'data/parcels/dallas_county_parcel_data/ACCOUNT_APPRL_YEAR.csv'
dallas_county_res_path = 'data/parcels/dallas_county_parcel_data/res_detail.csv'