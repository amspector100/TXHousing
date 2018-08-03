import time
import pandas as pd
import geopandas as gpd
from ..utilities import simple, spatial_joins
from . import boundaries

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
                                feature_columns_list = [houston_building_res_columns], county_level = False):
    """ Merge houston or harris parcel data with harris county data.

    :param feature_files: A list of paths of feature files to merge with the data. Defaults to
        the building_res file path. Can also be a string of one file (doesn't have to be a list). If you set this equal
        to None, it will just return the parcel shape data with no other data attached.
    :param feature_columns_list: A list of column headers for each file - you can find these in the
        Harris_Parcel_Feature_Columns microsoft database under data.
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
        place_shapes = place_shapes.loc[place_shapes['NAME'] == 'Houston']
        mapper = spatial_joins.fast_polygon_intersection(geodata, place_shapes, large_name_column='NAME')
        geodata['place'] = geodata.index.map(mapper)
        geodata = geodata.loc[geodata['place'].notnull()]

    # Return
    return geodata
