import time
import re # Regular expressions
import shapely
import numpy as np
import pandas as pd
import geopandas as gpd
import folium
from folium import LayerControl
import matplotlib.pyplot as plt
from plotnine import *

from inputs import *
from helpers import *
import spatial_functions as sf
import choropleth

# Global
time0 = time.time()

dallas_urban_areas = ['Dallas--Fort Worth--Arlington, TX', 'Denton--Lewisville, TX', 'McKinney, TX']
dallas_job_centers = ['Plano', 'Irving', 'Fort Worth', 'Arlington', 'Lewisville', 'McKinney', 'Rockwall', 'Garland', 'Denton', 'Frisco']

houston_urban_areas = ['Houston, TX', 'Conroe--The Woodlands, TX', 'Texas City, TX'] # I haven't run the processing function with the latter two as of 7/23/18
houston_job_centers = ['Rosenberg', 'Sugarland', 'The Woodlands', 'Katy', 'Pearland', 'La Porte', 'Friendswood']

austin_urban_areas = ['Austin, TX']
austin_job_centers = ['Round Rock', 'Georgetown', 'Cedar Park', 'Leandor', 'Taylor', 'Elgin', 'Bastrop', 'Lakeway']

def get_cached_parcel_path(name, level):
    return 'data/parcels/cached/{}_{}_parcels/parcels.shp'.format(name, level)

def get_cached_parcel_path_csv(name, level):
    return 'data/parcels/cached/{}_{}_parcels/parcels.csv'.format(name, level)

# Cache paths for Dallas
all_dallas_parcel_path = 'data/parcels/cached/dallas_all_parcels/parcels.shp'
dallas_metro_parcel_path = 'data/parcels/cached/dallas_metro_parcels/parcels.shp'
dallas_jobcenter_parcel_path = 'data/parcels/cached/dallas_jobcenter_parcels/parcels.shp'

# Cache paths for Houston
all_houston_parcel_path = 'data/parcels/cached/houston_all_parcels/parcels.shp'
houston_metro_parcel_path = 'data/parcels/cached/houston_metro_parcels/parcels.shp'
houston_jobcenter_parcel_path = 'data/parcels/cached/houston_jobcenter_parcels/parcels.shp'

# Cache paths for Austin
all_austin_parcel_path = 'data/parcels/cached/austin_all_parcels/parcels.shp'
austin_metro_parcel_path = 'data/parcels/cached/austin_metro_parcels/parcels.shp'
austin_jobcenter_parcel_path = 'data/parcels/cached/austin_jobcenter_parcels/parcels.shp'

# For analysis
step = 5
maximum = 45
outlier_maximum = 55

# Functions =------------------------------------------------------------------------------------------------

# Parse descriptions of lots and turn them into SF/MF/Other/Other Residential classifications
def process_lot_descriptions(gdf, description_column, broad_zone_dictionary):
    """
    :param gdf: A geodataframe, presumably of parcel data, but it could honestly be of anything.
    :param description_column: The column of descriptions to parse
    :param broad_zone_dictionary: A dictionary mapping base zones (i.e. 'Single Family') to a list of strings
    which signal that a description means that base zone. Ex: {'Single Family':['sf', 'single f'], 'Multifamily':['mf']}
    Note that the order of this dictionary DOES matter - the function will return the FIRST key which has a match.
    :return: A pandas series of the base zones.
    """
    # Create helper function
    def process_broad_zone(text):
        for key in broad_zone_dictionary:
            # This checks whether any of the strings in the broad_zone_dictionary[key] list appear in the text
            if bool(re.search(('|').join(broad_zone_dictionary[key]), text)):
                return key
        return "Other"

    broad_zones = gdf[description_column].apply(process_broad_zone)
    return broad_zones


def process_parcel_shapefile(path, zoning_input, county_name, broad_zone_dictionary, state_cd_feature, account_feature,
                             zipcode_feature = None, urban_areas_list = None, merge_paths = None,
                             left_keys = None, right_keys = None, gdf = None, fields_to_preserve = None, **kwargs):
    """
    :param path: Path of the shapefile to read in.
    :param zoning_input: The zoning input for the city of interest used to calculate distance from city center.
    :param county_name: Name of the county.
    :param broad_zone_dictionary: Dictionary mapping broad zones to keywords that we will use to convert the state_cd_feature
    into broad zones.
    :param state_cd_feature: The feature corresponding to state cd codes (or theoeretically could be another feature
    used to interpret broad zone, but state cds are preferred for consistency).
    :param account_feature: The feature corresponding to the tax account for each parcel.
    :param zipcode_feature: Default None. If the data already lists the zip codes of the parcels, then list the zipcode
    column name in this list to prevent unnecessary recalculation and spatial joins.
    :param merge_paths: If the geodata must be merged with another datasource, this should be an iterable containing the
    paths of the data to merge it with.
    :param left_keys: If the geodata must be merged with another datasource, this should be an iterable containing the
    left keys used to merge the data, in the same order of as the 'merge_paths'. Defaults to None, in which case the merge
    is performed using the account_feature as the merge key.
    :param right_keys: If the geodata must be merged with another datasource, this should be an iterable containing the
    right keys used to merge the data, in the same order of as the 'merge_paths'
    :param gdf: Optionally, supply a preread geodataframe instead of a path. Defaults to None.
    :param fields_to_preserve: A list of columns to preserve in the final data output.
    :param **kwargs: kwargs to pass to the read_csv call, if merging with external data.
    :return: A gdf, in lat long, with the following columns: geometry, area_sqft, broad_zone, county, place,
    ua, zipcode, account, state_cd, and any extra columns listed in the fields_to_preserve optional arg.
    """

    # Read in data
    if gdf is None:
        print('Reading {} parcel shapefile, time is {}'.format(county_name, time.time() - time0))
        gdf = gpd.read_file(path)

    # Transform to lat long, calculate area in square feet, and calculate distance to center
    print('Transforming {} parcel shapefile, time is {}'.format(county_name, time.time() - time0))
    gdf = gdf.loc[gdf['geometry'].apply(lambda x: x is not None)]
    gdf = gdf.loc[gdf['geometry'].is_valid]
    gdf = sf.get_area_in_units(gdf, scale = 1, name = 'area_sqft', final_projection = {'init':'epsg:4326'})
    gdf['dist_to_cent'] = sf.calculate_dist_to_center(gdf, zoning_input, drop_centroids = False)
    gdf['centroids'] = gdf['geometry'].centroid

    # Drop duplicate geometries by centroid
    gdf['centroids_string'] = gdf['centroids'].astype(str)
    gdf['long'] = gdf['centroids_string'].apply(lambda x: x.split('(')[1].split(' ')[0])
    gdf['lat'] = gdf['centroids_string'].apply(lambda x: x.split('(')[1].split(' ')[1].split(')')[0])
    gdf = gdf.drop_duplicates(subset=['centroids_string'], keep='first')


    # Merge with other data --
    if merge_paths is not None:
        print('Merging {} parcel shapefile, time is {}'.format(county_name, time.time() - time0))

        #  Start by ensuring that the left_keys, right_keys, merge_paths are iterables
        if isinstance(left_keys, str):
            left_keys = [left_keys]
        elif left_keys is None:
            left_keys = [account_feature]*len(merge_paths)
        if isinstance(right_keys, str):
            right_keys = [right_keys]
        if isinstance(merge_paths, str):
            merge_paths = [merge_paths]
        if len(merge_paths) != len(right_keys) or len(merge_paths) != len(left_keys):
            raise ValueError('The lengths of left_keys {}, right_keys {}, and merge_paths {} must be equal'.format(len(left_keys),
                                                                                                                    len(right_keys),
                                                                                                                    len(merge_paths)))

        # Read data and merge
        for merge_path, left_key, right_key in zip(merge_paths, left_keys, right_keys):
            extra_data = pd.read_csv(merge_path, **kwargs)
            extra_data = extra_data.loc[(extra_data[right_key].notnull())]
            extra_data = extra_data.drop_duplicates(subset = right_key, keep = 'first')
            gdf = gdf.merge(extra_data, left_on = left_key, right_on = right_key, how = 'left')

    # Parse base zone
    print('Processing nonspatial features in {} parcel shapefile, time is {}'.format(county_name, time.time() - time0))
    gdf[state_cd_feature] = gdf[state_cd_feature].replace([np.inf, -np.inf, np.nan], 'NONE').astype(str)
    broad_zone_dictionary['Unknown'] = ['NONE']
    gdf['broad_zone'] = process_lot_descriptions(gdf, state_cd_feature, broad_zone_dictionary)

    # Rename state_cd, account, and (potentially) zipcode fields
    try:
        gdf['account'] = gdf[account_feature].apply(lambda x: str(int(float(x))))
    except:
        gdf['account'] = gdf[account_feature].astype(str)

    gdf['state_cd'] = gdf[state_cd_feature].astype(str)
    if zipcode_feature is not None:
        gdf['zipcode'] = gdf[zipcode_feature].astype(str)

    # Now we start to do spatial calculations, using centroids for speed. Start by creating the spatial index and a
    # helper function.
    print('Processing spatial features in {} parcel shapefile, time is {}'.format(county_name, time.time() - time0))
    gdf['county'] = county_name
    gdf = gdf.set_geometry('centroids')
    spatial_index = gdf.sindex
    def Fast_Intersection(polygon, name, rowname, fragment_polygon = True, horiz = 10, vert = 10):
        """
        Checks which parcels intersect a polygon. This polygon presumably has a name (i.e. "Houston").
        Given a rowname (i.e. "place"), this function sets all of the parcels which intersect the polygon's rowname to
        the name of the polygon.
        :param polygon: The polygon of interest
        :param name: Name of hte polygon
        :param rowname: Rowname, i.e. place, to change in the gdf
        :param fragment_polygon: If true, fragment the polygon to speed up intersection.
        :return: None, but modifies the gdf. This function is capitalized because it has nonlocal (and therefore global)
        side effects intentionally.
        """
        nonlocal gdf
        print('For {} data, for {} geography, starting intersection for {}; time is {}'.format(county_name, rowname, name, time.time() - time0))
        if fragment_polygon:  # more than 15 square miles --> fragment
            polygon_list = sf.fragment(polygon, horiz = horiz, vert = vert)
        else:
            polygon_list = [polygon]

        all_precise_intersections = set()
        for grid_piece in polygon_list:
            possible_intersections_index = list(spatial_index.intersection(grid_piece.bounds))
            possible_intersections = gdf.iloc[possible_intersections_index]
            precise_intersections = set(possible_intersections.loc[possible_intersections.intersects(grid_piece)].index.tolist())
            all_precise_intersections = all_precise_intersections.union(precise_intersections)

        # Assign gdf the precise intersections indexes
        gdf.loc[all_precise_intersections, rowname] = name

    # Start by getting county polygon
    counties = gpd.read_file(county_boundaries_path)
    county_polygon = counties.loc[(counties['NAME'] == county_name) & (counties['STATEFP'] == '48'), 'geometry'].simplify(tolerance = 0.005).values[0]

    # Work through zip codes - begin by subsetting to the zip codes in the county
    if zipcode_feature is None:
        gdf['zipcode'] = "None"
        zipdata = get_zip_boundaries()
        zip_spatial_index = zipdata.sindex
        possible_intersections = zipdata.iloc[list(zip_spatial_index.intersection(county_polygon.bounds))]
        zipdata = possible_intersections.loc[possible_intersections.intersects(county_polygon).index.tolist()]

        # Now run the Fast_Intersection function
        for zipcode, zippolygon in zip(zipdata.index, zipdata['geometry']):
            Fast_Intersection(zippolygon, zipcode, 'zipcode')

    # Work thorough municipalities
    gdf['place'] = 'Unincorporated'
    places = gpd.read_file(texas_places_path)
    places_spatial_index = places.sindex
    possible_intersections = places.iloc[list(places_spatial_index.intersection(county_polygon.bounds))]
    places = places.loc[possible_intersections.intersects(county_polygon).index.tolist()]

    # Now run the Fast_Intersection function
    for placename, placepolygon in zip(places['NAME'], places['geometry']):
        Fast_Intersection(placepolygon, placename, 'place')

    gdf['ua'] = 'Other'
    if urban_areas_list is not None:
        uas = gpd.read_file(ua_path)
        uas = uas.loc[uas['NAME10'].isin(urban_areas_list)]
        for uaname, uapolygon in zip(uas['NAME10'], uas['geometry']):
            Fast_Intersection(uapolygon, uaname, 'ua', horiz = 20, vert = 20)

    # Subset and return
    final_column_list = ['area_sqft', 'broad_zone', 'county', 'place', 'ua', 'zipcode', 'account', 'state_cd', 'lat', 'long']
    if fields_to_preserve is not None:
        final_column_list.extend(fields_to_preserve)
    gdf = gpd.GeoDataFrame(gdf[final_column_list], geometry = gdf['geometry'])
    print(gdf)

    return gdf

# Take a big dataset (i.e. parcels, gdf) and subset to job centers
def subset_to_job_centers(gdf, job_centers, horiz = 3, vert = 3):
    """
    :param gdf: Some geodataframe
    :param job_centers: A list of names of job centers
    :return: The subset of the geodataframe which lies in the job centers, the job centers shapes,
    and the shapefile of places in texas.
    """

    # Get place shapes
    texas_places = gpd.read_file(texas_places_path)
    job_centers_shapes = texas_places.loc[texas_places['NAME'].isin(job_centers)]

    # Subset
    job_centers_dictionary = sf.fast_polygon_intersection(gdf, job_centers_shapes['geometry'], names = job_centers_shapes['NAME'], horiz = horiz, vert = vert)
    gdf['job_center'] = gdf.index.to_series().map(job_centers_dictionary)
    subsetted_data = gdf.loc[gdf['job_center'].notnull()]
    return subsetted_data, job_centers_shapes, texas_places

def subset_to_metro(gdf, names, geography = 'ua', places_to_ignore = None):
    """
    :param gdf: Some geodataframe.
    :param names: List of names to use to subset to a metro area
    :param geography: string, default 'ua'. Specifies the geography used to describe the metro - can either be 'ua',
    'csa', or 'cbsa'.
    :param places_to_ignore: List of places (places, not core based statistical areas) to ignore, i.e. ['Dallas']
    :return: The subset of the geodataframe lying in the 'names' csa/cbsa but not in the places_to_ignore.
    """

    # Get CSA/CBSA/UA areas
    if geography == 'csa':
        metro_shapes = gpd.read_file(csa_path)
    elif geography == 'cbsa':
        metro_shapes = gpd.read_file(cbsa_path)
    elif geography == 'ua':
        metro_shapes = gpd.read_file(ua_path)
        metro_shapes = metro_shapes.rename(columns = {'NAME10':'NAME'})
    else:
        raise ValueError("geography arg must be one of 'csa', 'cbsa', or 'ua', not {}".format(geography))

    boundaries = metro_shapes.loc[metro_shapes['NAME'].isin(names)]

    # Subset
    metro_dic = sf.fast_polygon_intersection(gdf, boundaries['geometry'], names = boundaries['NAME'], horiz =  25, vert = 25)
    gdf['metro'] = gdf.index.to_series().map(metro_dic)
    subsetted_data = gdf.loc[gdf['metro'].notnull()]

    # Get rid of core city to just analyze suburbs
    if places_to_ignore is not None:
        texas_places = gpd.read_file(texas_places_path)
        place_shapes = texas_places.loc[texas_places['NAME'].isin(places_to_ignore)]
        place_dic = sf.fast_polygon_intersection(subsetted_data, place_shapes['geometry'], names = place_shapes['NAME'], horiz = 10, vert = 10)
        subsetted_data['place'] = subsetted_data.index.to_series().map(place_dic).fillna(False)
        subsetted_data = subsetted_data.loc[subsetted_data['place'] == False]

    return subsetted_data, boundaries, metro_shapes

def get_all_dallas_parcel_data():

    tarrant_parcels = process_parcel_shapefile(path=processed_tarrant_county_parcel_path,
                                              zoning_input=dallas_inputs,
                                              county_name='Tarrant',
                                              broad_zone_dictionary=state_sptbcode_dictionary,
                                              state_cd_feature='Prprt_C',
                                              account_feature='TAXPIN',
                                              urban_areas_list=dallas_urban_areas,
                                              zipcode_feature=None,
                                              merge_paths=None,
                                              left_keys=None,
                                              right_keys=None,
                                              fields_to_preserve=None)

    # A3 refers to Condos, A4 refers to townhomes. Note I have not included property improvements (classes A6, A9, B6, B9)
    # - see the data dictionary in the data/parcels directory.
    collin_dictionary = {'Single Family':['A1'], 'Multifamily':['A3', 'A4', 'B1', 'B2', 'B3', 'B4']}
    collin_parcels = process_parcel_shapefile(path=collin_county_parcel_path,
                                              zoning_input=dallas_inputs,
                                              county_name='Collin',
                                              broad_zone_dictionary=collin_dictionary,
                                              state_cd_feature='state_cd',
                                              account_feature='PROP_ID',
                                              urban_areas_list=dallas_urban_areas,
                                              zipcode_feature=None,
                                              merge_paths=None,
                                              left_keys=None,
                                              right_keys=None,
                                              fields_to_preserve=None)

    denton_dictionary = {'Single Family': ['Single Family', 'Duplexes'],
                                   'Multifamily': ['Apartment', 'Townhomes', ', Condos'],
                                   'Other Residential': ['Residential']}
    denton_parcels = process_parcel_shapefile(path=denton_county_parcel_path,
                            zoning_input=dallas_inputs,
                            county_name='Denton',
                            broad_zone_dictionary=denton_dictionary,
                            state_cd_feature='CD_DESCRIP',
                            account_feature='PROP_ID',
                            urban_areas_list=dallas_urban_areas,
                            zipcode_feature=None,
                            merge_paths=None,
                            left_keys=None,
                            right_keys=None,
                            fields_to_preserve=None)

    dallas_parcels = process_parcel_shapefile(path=dallas_county_parcel_path,
                                              zoning_input=dallas_inputs,
                                              county_name='Dallas',
                                              broad_zone_dictionary=dallas_sptb_dictionary,
                                              state_cd_feature='SPTD_CD',
                                              account_feature='Acct',
                                              urban_areas_list=dallas_urban_areas,
                                              zipcode_feature=None,
                                              merge_paths=[dallas_county_land_path],
                                              left_keys=['Acct'],
                                              right_keys=['ACCOUNT_NUM'],
                                              fields_to_preserve=None)

    print('Combining and saving for Dallas, time is {}'.format(time.time() - time0))
    all_dallas_parcels = pd.concat([collin_parcels, dallas_parcels, denton_parcels, tarrant_parcels],
                                   axis=0,
                                   ignore_index=True)
    all_dallas_parcels.to_file(all_dallas_parcel_path)
    all_dallas_parcels[[col for col in all_dallas_parcels.columns if col != 'geometry']].to_csv(get_cached_parcel_path_csv('dallas', 'all'))
    print('Finished with Dallas, time is {}'.format(time.time() - time0))




def get_all_austin_parcel_data():

    # See https://tax-office.traviscountytx.gov/pages/SPTC.php - A4/A5 refer to condos.
    travis_dictionary = {'Single Family':['A1', 'A2', 'A3'], 'Multifamily':['A4', 'A5', 'B1', 'B2', 'B3', 'B4']}
    travis_parcels = process_parcel_shapefile(path=travis_county_parcel_path,
                                                 zoning_input=austin_inputs,
                                                 county_name='Travis',
                                                 broad_zone_dictionary=travis_dictionary,
                                                 state_cd_feature='state_cd',
                                                 account_feature='PROP_ID',
                                                 urban_areas_list=austin_urban_areas,
                                                 zipcode_feature=None,
                                                 merge_paths=[travis_county_data_path],
                                                 left_keys=['PROP_ID'],
                                                 right_keys=['PROP_ID'],
                                                 fields_to_preserve=None)

    # See https://www.wcad.org/wp-content/uploads/2016/09/2015Report.pdf for the dictionary information
    williamson_dictionary = {'Single Family':['A1', 'A9'], 'Multifamily':['A8', 'B1', 'B2', 'B4']}
    williamson_parcels = process_parcel_shapefile(path=williamson_county_parcel_path,
                                                 zoning_input=austin_inputs,
                                                 county_name='Williamson',
                                                 broad_zone_dictionary=williamson_dictionary,
                                                 state_cd_feature='StateCode',
                                                 account_feature='PIN',
                                                 urban_areas_list=austin_urban_areas,
                                                 zipcode_feature=None,
                                                 merge_paths=[williamson_county_real_improvement_path],
                                                 left_keys=['PIN'],
                                                 right_keys=['QuickRefID'],
                                                 fields_to_preserve=None)

    print('Combining and saving for Austin, time is {}'.format(time.time() - time0))
    all_austin_parcels = pd.concat([travis_parcels, williamson_parcels], axis=0, ignore_index=True)
    all_austin_parcels = all_austin_parcels.loc[~all_austin_parcels['geometry'].apply(lambda x: x is None)]
    all_austin_parcels = all_austin_parcels.loc[all_austin_parcels['geometry'].is_valid]
    all_austin_parcels.to_file(all_austin_parcel_path) # shapefile
    all_austin_parcels[[col for col in all_austin_parcels.columns if col != 'geometry']].to_csv(get_cached_parcel_path_csv('austin', 'all')) # csv
    print('Finished with all Austin parcels, time is {}'.format(time.time() - time0))

def get_all_houston_parcel_data():

    # Harris County -- here, 1006 refers to condominiums, 1007 refers to townhomes
    #harris_dictionary = {'Single Family': ['1001'], 'Multifamily': ['1002', '1003', '1004', '1005', '1006', '1007', '4209', '4211', '4212', '4214', '4299']}
    harris_parcels = process_houston_parcel_data(county_level = True)
    harris_parcels = process_parcel_shapefile(path = None,
                             zoning_input = houston_inputs,
                             county_name = 'Harris',
                             broad_zone_dictionary = state_sptbcode_dictionary,
                             state_cd_feature = 'USE_CODE',
                             account_feature = 'HCAD_NUM',
                             urban_areas_list=houston_urban_areas,
                             zipcode_feature=None,
                             merge_paths=None,
                             left_keys=None,
                             right_keys=None,
                             fields_to_preserve=None,
                             gdf = harris_parcels)

    # Fort Bend county
    fort_bend_parcels = process_parcel_shapefile(path = fort_bend_parcel_path,
                             zoning_input = houston_inputs,
                             county_name = 'Fort Bend',
                             broad_zone_dictionary = state_sptbcode_dictionary,
                             state_cd_feature = 'LMainSegSP',
                             account_feature = 'UID',
                             urban_areas_list=houston_urban_areas,
                             zipcode_feature=None,
                             merge_paths=None,
                             left_keys=None,
                             right_keys=None,
                             fields_to_preserve=None)

    # Montgomery county
    montgomery_parcels= process_parcel_shapefile(path = montgomery_county_parcel_path,
                             zoning_input = houston_inputs,
                             county_name = 'Montgomery',
                             broad_zone_dictionary = state_sptbcode_dictionary,
                             state_cd_feature = 'fStateCode',
                             account_feature = 'PropertyNu',
                             urban_areas_list = houston_urban_areas,
                             zipcode_feature=None,
                             merge_paths=None,
                             left_keys=None,
                             right_keys=None,
                             fields_to_preserve=None)


    print('Combining and saving at {}'.format(time.time() - time0))
    all_houston_parcels = pd.concat([harris_parcels, fort_bend_parcels, montgomery_parcels], axis=0, ignore_index=True)
    all_houston_parcels.to_file(all_houston_parcel_path)
    all_houston_parcels[[col for col in all_houston_parcels.columns if col != 'geometry']].to_csv(get_cached_parcel_path_csv('houston', 'all')) # csv

    print('Finished saving, time is {}'.format(time.time() - time0))

# Helper functions which generate outfile paths for lotsize/zoning data
def get_lotsize_path(name, level):
    return 'data/caches/suburbs/{}_{}_lotsizes.csv'.format(name, level)

def get_percent_zoned_path(name, level):
    return 'data/caches/suburbs/{}_{}_parcel_zone_rings.csv'.format(name, level)

def get_municipality_choropleth_path(name):
    return 'data/caches/suburbs/{}_municipality_calculations/municipality_calculations.shp'.format(name)

# Commuting ------------------------------------------------------------------------------------------------------
def analyze_transportation_networks(names, zoning_inputs, num_blocks = 5000, step = 2.5, maximum = 60):

    block_data = sf.get_block_geodata(['X08_COMMUTING', 'X01_AGE_AND_SEX'], cities = None)

    # Get total/local workers
    block_data['local_workers'] = block_data['B08008e3'] + block_data['B08008e8'] # Female and male workers working in their place of residence
    block_data['total_workers'] = block_data['B08008e2'] + block_data['B08008e7'] # Total number of male and female workers in the block group


    # Create spatial indexes and subset
    spatial_index = block_data.sindex
    all_data = gpd.GeoDataFrame()
    for name, zoning_input in zip(names, zoning_inputs):

        # Query and find nearest neighbors, subset
        nearest_index = list(spatial_index.nearest((zoning_input.long, zoning_input.lat), num_results=num_blocks))
        city_data = block_data.iloc[nearest_index]
        city_data['dist_to_center'] = sf.calculate_dist_to_center(city_data, zoning_input)
        city_data['City'] = name
        all_data = pd.concat([all_data, city_data])

    # Average commute time

    # Old way of calculating --- my way is better ---
    #all_data = all_data.loc[all_data['B08135e1'] != 0]
    # all_data['avg_commute_time'] = all_data['B08135e1'].divide(all_data['B08303e1'])

    all_data['smoothed_dist_to_center'] =  all_data['dist_to_center'].apply(lambda x: step*(np.round(x/step)))
    commute_brackets = {'B08303e2':2.5,
                        'B08303e3':7.5,
                        'B08303e4':12.5,
                        'B08303e5':17.5,
                        'B08303e6':22.5,
                        'B08303e7':27.5,
                        'B08303e8':32.5,
                        'B08303e9':37.5,
                        'B08303e10':42.5,
                        'B08303e11':52.5,
                        'B08303e12':75,
                        'B08303e13':90}
    all_data['total_commute_time'] = 0
    for key in commute_brackets:
        all_data['total_commute_time'] += all_data[key]*commute_brackets[key]
    print(all_data[['total_commute_time', 'B08135e1']])

    # Calculate averages in rings
    total_commuters = all_data.groupby(['smoothed_dist_to_center', 'City'])['B08303e1'].sum()
    total_commute_time = all_data.groupby(['smoothed_dist_to_center', 'City'])['total_commute_time'].sum()
    result = pd.DataFrame(total_commute_time.divide(total_commuters))
    result.reset_index(inplace = True)
    result = result.rename(columns = {0:'avg_commute_time'})
    result = result.loc[result['smoothed_dist_to_center'] <= maximum]

    plot = (ggplot(result, aes(x = 'smoothed_dist_to_center', y = 'avg_commute_time', fill = 'City'))
                   + geom_col(position = 'dodge')
                   + facet_wrap('~City')
                   + theme_bw()
                   + labs(title = 'Commute Times by Distance from the City Center in Austin, Dallas, and Houston',
                          x = 'Distance from City Center (Miles)', y = 'Average Commute Time (Minutes)'))
    plot.save('Figures/Suburbs/travel_times.svg', width = 12, height = 10)

    total_workers = all_data.groupby(['smoothed_dist_to_center', 'City'])['total_workers'].sum()
    total_local_workers = all_data.groupby(['smoothed_dist_to_center', 'City'])['local_workers'].sum()
    workers_pct = 100*pd.DataFrame(total_local_workers.divide(total_workers))
    workers_pct.reset_index(inplace = True)
    workers_pct = workers_pct.rename(columns = {0:'local_workers_pct'})
    workers_pct = workers_pct.loc[workers_pct['smoothed_dist_to_center'] <= maximum]

    plot = (ggplot(workers_pct, aes(x = 'smoothed_dist_to_center', y = 'local_workers_pct', fill = 'City'))
                   + geom_col(position = 'dodge')
                   + facet_wrap('~City')
                   + theme_bw()
                   + labs(title = 'Percent of Workers Working in Place of Residence by Distance from City Center',
                          x = 'Distance from City Center (Miles)', y = 'Percent of Workers Working in Place of Residence'))
    plot.save('Figures/Suburbs/local_workers.svg', width = 12, height = 10)


    #city_data['avg_commute_time'] = city_data['B08135e1'].divide(city_data['B08303e1'])



def analyze_dallas_suburbs(job_centers = dallas_job_centers):
    time0 = time.time()

    # Analyze the big north texas zoning shapefile --------------------------------------

    # Basic analysis of percent zoned SF in the suburbs
    north_texas_data = process_zoning_shapefile(north_texas_inputs)
    ua_zones, boundaries, ua_shapes = subset_to_metro(north_texas_data, dallas_urban_areas, geography = 'ua', places_to_ignore = ['Dallas'])

    print('Starting to do the polygons_intersect_rings')
    ua_result = sf.polygons_intersect_rings(ua_zones, dallas_inputs, factor = 'broad_zone', step = step, maximum = maximum, outlier_maximum = outlier_maximum, categorical = True)
    ua_result.to_csv('data/caches/suburbs/dallas_ua_zone_rings.csv')
    print('Time is {}'.format(time.time() - time0))

    job_center_zones, job_centers_shapes, texas_places = subset_to_job_centers(north_texas_data, job_centers)
    job_result = sf.polygons_intersect_rings(job_center_zones, dallas_inputs, factor = 'broad_zone', step = step, maximum = maximum, outlier_maximum = outlier_maximum, categorical = True)
    job_result.to_csv('data/caches/suburbs/dallas_jobcenter_zone_rings.csv')
    print('Time is {}'.format(time.time() - time0))

def plot_municipality_choropleth(name, zoning_input, job_centers, to_exclude = None):

    path = get_municipality_choropleth_path(name)
    basemap = folium.Map([zoning_input.lat, zoning_input.long], zoom_start = 9)

    # Read and format data
    geodata = gpd.read_file(path)
    geodata = geodata.set_index('NAME')
    geodata["NAME"] = geodata.index
    geodata = geodata.loc[(geodata['Single F_1'] != 0) & (geodata['Other-zone'] != 0)]
    if to_exclude is not None:
        geodata = geodata.loc[[not bool for bool in geodata['NAME'].isin(to_exclude)]]

    geodata = geodata.rename(columns = {'Other-zone':'Percent of Land Developed as Nonresidential',
                                        'Other-lots':'Other-lotsize',
                                        'Multifamil':'Multifamily-lotsize',
                                        'Multifam_1':'Percent of Land Used as Multifamily',
                                        'Single Fam':'Single Family Average Lotsize (Square Feet)',
                                        'Single F_1':'Percent of Land Used as Single Family'})
    geodata.crs = {'init':'epsg:4326'}

    # Get rid of places for which no parcels intersected
    for factor in ['Percent of Land Used as Single Family', 'Percent of Land Used as Multifamily', 'Percent of Land Developed as Nonresidential', 'Single Family Average Lotsize (Square Feet)']:
        geojson, colormap = choropleth.continuous_choropleth(geodata, factor, factor, scale_name = None, colors = ['lightgreen', 'navy', 'black'],
                                                             show = False, geometry_column = 'geometry', basemap = basemap)

    for jc in job_centers:
        try:
            coords = choropleth.retrieve_coords(geodata.loc[geodata["NAME"] == jc, 'geometry'].values[0].centroid)
            folium.Marker(
                location=coords,
                popup=jc,  # The name
                icon=folium.Icon(color='gray')
            ).add_to(basemap)
        except IndexError: # Occurs if the jobcenter is no longer in the place shape, presumably because we had no parcel data on it
            continue


    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    folium.TileLayer('CartoDB positron').add_to(basemap)
    folium.TileLayer('Stamen Toner').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save('Figures/Suburbs/{}_suburb_choropleth.html'.format(name))

def analyze_land_use_by_metro(name):
    print(name)

    # Read data
    def land_use_by_municipality_path(name):
        return 'data/caches/suburbs/{}_land_use_by_municipality.csv'.format(name)

    path = get_cached_parcel_path_csv(name, 'all')
    data = pd.read_csv(path, engine = 'python')
    data = data.drop_duplicates(subset = ['lat', 'long'], keep = 'first')

    # Calculate percent of area zoned
    zone_areas = data.groupby(['broad_zone', 'place'])['area_sqft'].sum()
    municipality_areas = data.groupby(['place'])['area_sqft'].sum()
    final_data = zone_areas.divide(municipality_areas)
    final_data = final_data.reset_index()
    final_data.columns = ['broad_zone', 'place', 'Percent of Land Used']
    final_data.to_csv(land_use_by_municipality_path(name))

def plot_choropleth_close_to_point(lat, long, zoning_input, gdf, spatial_index, save_path, column = None, num_polygons = 2500):

    basemap = folium.Map([zoning_input.lat, zoning_input.long], zoom_start = 9)
    ids = list(spatial_index.nearest([long, lat], num_results = num_polygons))
    subset = gdf.iloc[ids]
    subset.crs = {'init':'epsg:4326'}
    if column is not None:
        choropleth.categorical_choropleth(subset, column).add_to(basemap)
    else:
        subset.loc[:, 'dummy'] = 0
        choropleth.categorical_choropleth(subset, 'dummy').add_to(basemap)
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    folium.TileLayer('CartoDB positron').add_to(basemap)
    basemap.save(save_path)




if __name__ == '__main__':
    #get_all_houston_parcel_data()

    #analyze_transportation_networks(['austin', 'dallas', 'houston'], [austin_inputs, dallas_inputs, houston_inputs])
    #plot_municipality_choropleth('austin', austin_inputs, austin_job_centers)
    #plot_municipality_choropleth('dallas', dallas_inputs, dallas_job_centers, to_exclude = ['Combine'])
    #plot_municipality_choropleth('houston', houston_inputs, houston_job_centers)

    for cityname, cityinput, universitylats, universitylongs, universitynames in zip(['austin', 'dallas', 'houston'],
                                                                                     [austin_inputs, dallas_inputs, houston_inputs],
                                                                                     [[39.2846914], [32.8412, 32.7299], [29.7174, 29.7199]],
                                                                                     [[-97.7416248], [-96.7845, -97.1140], [-95.4018, -95.3422]],
                                                                                     [['UT_Austin'], ['Southern_Methodist', 'UT_Arlington'], ['Rice', 'UHouston']]):

        print(cityname)
        ou_savepath = 'Figures/Testing/{}_core_others_and_unknown_parcels.html'.format(cityname)
        data = gpd.read_file(get_cached_parcel_path(cityname, 'all'))
        print('Finished reading at time {}'.format(time.time() - time0))
        spatind = data.sindex
        print('Found sindex at time {}'.format(time.time() - time0))
        others_and_unknowns = data.loc[data['broad_zone'].isin(['Other', 'Unknown'])]
        plot_choropleth_close_to_point(cityinput.lat, cityinput.long, cityinput, others_and_unknowns, others_and_unknowns.sindex,
                                       save_path = ou_savepath, num_polygons = 4000)

        for lat, long, uname in zip(universitylats, universitylongs, universitynames):
            print(uname)
            uni_save_path = 'Figures/Testing/{}_surrounding_parcels.html'.format(uname)
            plot_choropleth_close_to_point(lat, long, cityinput, data, spatind, save_path = uni_save_path, column = 'broad_zone', num_polygons = 2500)

