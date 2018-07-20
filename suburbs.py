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
state_sptbcode_dictionary = {'Single Family':['A1', 'A2'], 'Multifamily':['A3', 'A4', 'B1', 'B2', 'B3', 'B4']} # Used for a lot of different counties

dallas_urban_areas = ['Dallas--Fort Worth--Arlington, TX', 'Denton--Lewisville, TX', 'McKinney, TX']
dallas_job_centers = ['Plano', 'Irving', 'Fort Worth', 'Arlington', 'Lewisville', 'McKinney', 'Rockwall', 'Garland', 'Denton', 'Frisco']

houston_urban_areas = ['Houston, TX']
houston_job_centers = ['Rosenberg', 'Sugarland', 'The Woodlands', 'Katy', 'Pearland', 'La Porte', 'Friendswood']

austin_urban_areas = ['Austin, TX']
austin_job_centers = ['Round Rock', 'Georgetown', 'Cedar Park', 'Leandor', 'Taylor', 'Elgin', 'Bastrop', 'Lakeway']

def get_cached_parcel_path(name, level):
    return 'data/parcels/cached/{}_{}_parcels/parcels.shp'.format(name, level)

def get_cached_parcel_path_csv(name, level):
    return 'data/parcels/cached/{}_{}_parcels/parcels.csv'.format(name, level)

# Cache paths for Dallas
all_dallas_parcel_path = 'data/parcels/cached/all_dallas_region_parcels/parcels.shp'
dallas_metro_parcel_path = 'data/parcels/cached/dallas_metro_parcels/parcels.shp'
dallas_jobcenter_parcel_path = 'data/parcels/cached/dallas_jobcenter_parcels/parcels.shp'

# Cache paths for Houston
all_houston_parcel_path = 'data/parcels/cached/all_houston_region_parcels/parcels.shp'
houston_metro_parcel_path = 'data/parcels/cached/houston_metro_parcels/parcels.shp'
houston_jobcenter_parcel_path = 'data/parcels/cached/houston_jobcenter_parcels/parcels.shp'

# Cache paths for Austin
all_austin_parcel_path = 'data/parcels/cached/all_austin_region_parcels/parcels.shp'
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

def process_parcel_shapefile(path, feature, broad_zone_dictionary, name = '', dropnas = False, gdf = None):

    # Read data unless data is supplied
    if gdf is None:
        print('Reading {} parcel shapefile, time is {}'.format(name, time.time() - time0))
        gdf = gpd.read_file(path)

    # Drop nas if necessary
    if dropnas:
        gdf[feature] = gdf[feature].replace([np.inf, -np.inf], np.nan)
        gdf = gdf.dropna(subset=[feature], how='all')


    # Get basezone
    print('Processing {} parcel data, time is {}'.format(name, time.time() - time0))
    gdf[feature] = gdf[feature].astype(str)
    gdf['broad_zone'] = process_lot_descriptions(gdf, feature, broad_zone_dictionary)
    gdf.drop([col for col in gdf.columns if col not in ['broad_zone', 'geometry']], axis = 1, inplace = True)

    # Process to get rid of empty and invalid geometries
    gdf = gdf.loc[[not bool for bool in gdf['geometry'].apply(lambda x: x is None)]]
    gdf = gdf.loc[gdf['geometry'].is_valid]

    # Transform
    print('Transforming {} parcel data, time is {}'.format(name, time.time() - time0))
    try:
        gdf = gdf.to_crs({'init':'epsg:4326'})
    except Exception as e:
        print(e)
        print('Failed to transform to lat/long, here is the name, crs and geometry')
        print(gdf['geometry'])
        print(name, gdf.crs)
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

def get_all_dallas_parcel_data(combine_parcel_data = False, subset_parcel_data = False):


    time0 = time.time()

    if combine_parcel_data:

        # Do dallas manually because it requires a join with non spatial data
        # Dallas --
        print('Reading Dallas county parcels, time is {}'.format(time.time() - time0))
        dallas_county_parcels = gpd.read_file(dallas_county_parcel_path) # 2018
        print('Processing Dallas county parcels broad zones, time is {}'.format(time.time() - time0))
        dallas_land_data = pd.read_csv(dallas_county_land_path) # 2018
        dallas_land_data = dallas_land_data.rename(columns = {'ACCOUNT_NUM':'Acct'}).drop_duplicates(subset = 'Acct', keep = 'first')
        dallas_land_data['ZONING'] = dallas_land_data['ZONING'].astype(str)
        dallas_county_parcels = dallas_county_parcels.merge(dallas_land_data, on = 'Acct')
        dallas_base_zone_dictionary = {'Single Family':['SINGLE F', 'SF', 'DUPLEX', 'SINGLE-F', 'ONE FAM'],
                                       'Multifamily':['MIXED USE', 'MULTIF', 'MULTI-F', 'MULTI F', 'MULTIPLE-F',
                                                      'MULTIPLE F', 'APARTM', 'MF-', 'TOWNH', 'TOWN H', 'QUADRAPLEX',
                                                      'FOUR HOME', 'TWO F'],
                                       'Other Residential':['RESIDENTIAL', 'DWELLING']}
        dallas_county_parcels['broad_zone'] = process_lot_descriptions(dallas_county_parcels, 'ZONING', dallas_base_zone_dictionary)
        dallas_county_parcels.drop([col for col in dallas_county_parcels.columns if col not in ['broad_zone', 'geometry']], inplace = True, axis = 1)  # Save memory by subsetting
        dallas_county_parcels = dallas_county_parcels.to_crs({'init':'epsg:4326'})

        # Denton --
        denton_base_zone_dictionary = {'Single Family':['Single Family', 'Duplexes'], 'Multifamily':['Apartment', 'Townhomes', ', Condos'],
                                       'Other Residential':['Residential']}
        denton_county_parcels = process_parcel_shapefile(denton_county_parcel_path, 'CD_DESCRIP', denton_base_zone_dictionary, name = 'Denton')

        # Collin --
        # A3 refers to Condos, A4 refers to townhomes. Note I have not included property improvements (classes A6, A9, B6, B9)
        # - see the data dictionary in the data/parcels directory.
        collin_base_zone_dictionary = state_sptbcode_dictionary
        collin_county_parcels = process_parcel_shapefile(collin_county_parcel_path, 'state_cd', collin_base_zone_dictionary, name = 'Collin')
        print(collin_county_parcels)

        # Tarrant --
        tarrant_base_zone_dictionary = collin_base_zone_dictionary # Same state codes
        tarrant_county_parcels = process_parcel_shapefile(processed_tarrant_county_parcel_path, 'Prprt_C', tarrant_base_zone_dictionary, name = 'Tarrant')
        print(tarrant_county_parcels)

        # Combine and cache
        all_dallas_parcels = pd.concat([denton_county_parcels, dallas_county_parcels, tarrant_county_parcels, collin_county_parcels], axis = 0, ignore_index = True)
        all_dallas_parcels.to_file(all_dallas_parcel_path)

    # Geographically subset
    if subset_parcel_data:

        print('Reading cached (combined) parcel data, time is {}'.format(time.time() - time0))
        all_dallas_parcels = gpd.read_file(all_dallas_parcel_path)

        print('Subsetting to metro, time is {}'.format(time.time() - time0))
        metro_dallas_parcels, metro_shapes, texas_uas = subset_to_metro(all_dallas_parcels, dallas_urban_areas, geography = 'ua', places_to_ignore = ['Dallas'])
        metro_dallas_parcels.drop('centroids', inplace = True, axis = 1)
        metro_dallas_parcels.to_file(dallas_metro_parcel_path)

        print('Subsetting to jobcenters, time is {}'.format(time.time() - time0))
        jobcenter_dallas_parcels, jobcenter_shapes, texas_places = subset_to_job_centers(all_dallas_parcels, dallas_job_centers)
        jobcenter_dallas_parcels.drop('centroids', inplace = True, axis = 1)
        jobcenter_dallas_parcels.to_file(dallas_jobcenter_parcel_path)

    print('Finished with Dallas parcels, time is {}'.format(time.time() - time0))


def get_all_austin_parcel_data(combine_parcel_data = False, subset_parcel_data = True):

    if combine_parcel_data:

        # Williamson -- already in lat/long
        print('Reading Williamson parcels, time is {}'.format(time0 - time.time()))
        williamson_parcels = gpd.read_file(williamson_county_parcel_path)
        williamson_data = pd.read_csv(williamson_county_real_improvement_path)
        williamson_data.drop_duplicates(subset='PropertyID', inplace=True)
        williamson_parcels = williamson_parcels.merge(williamson_data, on = 'PropertyID', how = 'left')
        williamson_parcels['fSPTB'] = williamson_parcels['fSPTB'].astype(str)
        williamson_dictionary = state_sptbcode_dictionary
        williamson_parcels['broad_zone'] = process_lot_descriptions(williamson_parcels, 'fSPTB', williamson_dictionary)
        print(williamson_parcels[['broad_zone', 'fSPTB']])
        williamson_parcels.drop([col for col in williamson_parcels.columns if col not in ['broad_zone', 'geometry']], inplace = True, axis = 1)

        # Travis --
        print('Reading Travis parcels, time is {}'.format(time.time() - time0))
        travis_parcels = gpd.read_file(travis_county_parcel_path)
        travis_data = pd.read_csv(travis_county_data_path)
        travis_data.drop_duplicates(subset = 'PROP_ID', inplace = True)
        travis_parcels = travis_parcels.merge(travis_data, on = 'PROP_ID', how = 'left')
        travis_parcels['state_cd'] = travis_parcels['state_cd'].astype(str)
        travis_dictionary = state_sptbcode_dictionary
        travis_parcels['broad_zone'] = process_lot_descriptions(travis_parcels, 'state_cd', travis_dictionary)
        print('Percent of state_cds which are not null is {}'.format(travis_parcels['state_cd'].notnull().sum()/travis_parcels.shape[0]))
        print(travis_parcels[['broad_zone', 'state_cd']])
        travis_parcels.drop([col for col in travis_parcels.columns if col not in ['broad_zone', 'geometry']], inplace = True, axis = 1)
        travis_parcels = travis_parcels.loc[[not bool for bool in travis_parcels['geometry'].apply(lambda x: x is None)]]
        travis_parcels = travis_parcels.loc[travis_parcels['geometry'].is_valid]
        travis_parcels = travis_parcels.to_crs({'init':'epsg:4326'})

        all_austin_parcels = pd.concat([travis_parcels, williamson_parcels], axis=0, ignore_index=True)
        all_austin_parcels = all_austin_parcels.loc[[not bool for bool in all_austin_parcels['geometry'].apply(lambda x: x is None)]]
        all_austin_parcels = all_austin_parcels.loc[all_austin_parcels['geometry'].is_valid]

        print(all_austin_parcels)
        all_austin_parcels.to_file(all_austin_parcel_path)

    if subset_parcel_data:

        all_austin_parcels = gpd.read_file(all_austin_parcel_path)
        print(all_austin_parcels)

        print('Subsetting to metro, time is {}'.format(time.time() - time0))
        austin_metro_parcels, metro_shapes, texas_uas = subset_to_metro(all_austin_parcels, austin_urban_areas, geography = 'ua', places_to_ignore = ['Austin'])
        austin_metro_parcels.drop('centroids', inplace = True, axis = 1)
        austin_metro_parcels.to_file(austin_metro_parcel_path)

        print('Subsetting to jobcenters, time is {}'.format(time.time() - time0))
        austin_jobcenter_parcels, jobcenter_shapes, texas_places = subset_to_job_centers(all_austin_parcels, austin_job_centers)
        austin_jobcenter_parcels.drop('centroids', inplace = True, axis = 1)
        austin_jobcenter_parcels.to_file(austin_jobcenter_parcel_path)

        print('Finished with Austin parcels, time is {}'.format(time.time() - time0))


    print('Finished, time is {}'.format(time.time() - time0))


def get_all_houston_parcel_data(combine_parcel_data = False, subset_parcel_data = False):


    # Combine and process parcel data
    if combine_parcel_data:

        # Montgomery county
        montgomery_dictionary = state_sptbcode_dictionary
        montgomery_parcels = process_parcel_shapefile(montgomery_county_parcel_path, 'fStateCode', montgomery_dictionary, name = 'Montgomery')
        print(montgomery_parcels)

        # Fort Bend county
        fort_bend_dictionary = state_sptbcode_dictionary
        fort_bend_parcels = process_parcel_shapefile(fort_bend_parcel_path, 'LMainSegSP', fort_bend_dictionary, name = 'Fort Bend')
        print(fort_bend_parcels)

        # Harris County -- here, 1006 refers to condominiums, 1007 refers to townhomes
        harris_dictionary = {'Single Family': ['1001'], 'Multifamily': ['1002', '1003', '1004', '1005', '1006', '1007', '4209', '4211', '4212', '4214', '4299']}
        harris_parcels = process_houston_parcel_data(county_level = True)
        harris_parcels = harris_parcels.loc[harris_parcels['LAND_USE_CODE'].notnull()]
        harris_parcels['LAND_USE_CODE'] = harris_parcels['LAND_USE_CODE'].apply(lambda x: str(int(x)))
        harris_parcels = process_parcel_shapefile(None, 'LAND_USE_CODE', harris_dictionary, name = 'Harris', dropnas = True, gdf = harris_parcels)
        print(harris_parcels)

        all_houston_parcels = pd.concat([harris_parcels, fort_bend_parcels, montgomery_parcels], axis=0, ignore_index=True)
        print(all_houston_parcels)
        all_houston_parcels.to_file(all_houston_parcel_path)

    if subset_parcel_data:

        all_houston_parcels = gpd.read_file(all_houston_parcel_path)
        print(all_houston_parcels)

        print('Subsetting to metro, time is {}'.format(time.time() - time0))
        houston_metro_parcels, metro_shapes, texas_uas = subset_to_metro(all_houston_parcels, houston_urban_areas, geography = 'ua', places_to_ignore = ['Houston'])
        houston_metro_parcels.drop('centroids', inplace = True, axis = 1)
        houston_metro_parcels.to_file(houston_metro_parcel_path)

        print('Subsetting to jobcenters, time is {}'.format(time.time() - time0))
        jobcenter_houston_parcels, jobcenter_shapes, texas_places = subset_to_job_centers(all_houston_parcels, houston_job_centers)
        jobcenter_houston_parcels.drop('centroids', inplace = True, axis = 1)
        jobcenter_houston_parcels.to_file(houston_jobcenter_parcel_path)

        print('Finished with Houston parcels, time is {}'.format(time.time() - time0))

# Helper functions which generate outfile paths for lotsize/zoning data
def get_lotsize_path(name, level):
    return 'data/caches/suburbs/{}_{}_lotsizes.csv'.format(name, level)

def get_percent_zoned_path(name, level):
    return 'data/caches/suburbs/{}_{}_parcel_zone_rings.csv'.format(name, level)

# Analyze the % of land zoned as SF as well as the lot sizes of the parcel data
def analyze_suburb_parcel_data(metro_path, jobcenters_path, name, zoning_input):
    """
    :param general_path: The path of the processed parcels which lie in the urban area surrounding the municipality but do not lie inside the municipality itself.
    :param jobcenters_path: The path of the processed parcels which lie in the job centers.
    :param name: Name of the city (used for creating outfile paths)
    :param zoning_input: The zoning_input for the city. Used to find the city center.
    :return: None. More importantly, this writes the results to csvs under data/caches/suburbs.
    """

    # Helper function
    def analyze_suburb_data_subset(gdf, level):

        # Get percent zoned - this will take a while --
        percent_zoned = sf.polygons_intersect_rings(gdf,
                                                    zoning_input,
                                                    factor = 'broad_zone',
                                                    categorical=True,
                                                    step = step, maximum = maximum,
                                                    outlier_maximum = outlier_maximum)

        percent_zoned_path = get_percent_zoned_path(name, level)
        percent_zoned.to_csv(percent_zoned_path)

        # Get lotsizes - this will also take a while --
        # Get area of parcels in square feet - this is expensive and will take a while
        gdf = sf.get_area_in_units(gdf, scale = 1, name = 'area')

        # For efficiency, create centroids and use these
        gdf['centroids'] = gdf['geometry'].centroid
        gdf = gdf.set_geometry('centroids')

        # Subset and initialize result (subsetting saves time)
        gdf = gdf.loc[gdf['broad_zone'].isin(['Single Family', 'Multifamily'])]
        lotsizes = pd.DataFrame()

        for zone in ['Single Family', 'Multifamily']:
            lotsizes[zone] = sf.points_intersect_rings(gdf.loc[gdf['broad_zone'] == zone],
                                                     zoning_input,
                                                     factor='area',
                                                     categorical=False,
                                                     by='mean',
                                                     per_square_mile=False,
                                                     geometry_column='centroids',
                                                     step = step,
                                                     maximum = maximum)

        lotsizes_path = get_lotsize_path(name, level)
        lotsizes.to_csv(lotsizes_path)

    print('Reading jobs data for {} at time {}'.format(name, time.time() - time0))
    jobs_data = gpd.read_file(jobcenters_path)
    print('Analyzing jobs data for {} at time {}'.format(name, time.time() - time0))
    analyze_suburb_data_subset(jobs_data, 'jobcenter')
    del jobs_data
    print('Reading ua data for {} at time {}'.format(name, time.time() - time0))
    ua_data = gpd.read_file(metro_path)
    print('Analyzing ua data for {} at time {}'.format(name, time.time() - time0))
    analyze_suburb_data_subset(ua_data, 'ua')
    del ua_data
    print('Finished for {} at time {}'.format(name, time.time() - time0))
    return None

# Just runs the analysis for everything
def analyze_all_suburb_parcel_data():

    for metro_path, jobcenters_path, name, zoning_input in zip([austin_metro_parcel_path,
                                                                dallas_metro_parcel_path,
                                                                houston_metro_parcel_path],
                                                               [austin_jobcenter_parcel_path,
                                                                dallas_jobcenter_parcel_path,
                                                                houston_jobcenter_parcel_path],
                                                               ['Austin', 'Dallas', 'Houston'],
                                                               [austin_inputs, dallas_inputs, houston_inputs]):

        print('Starting on {} at time {}'.format(name, time.time() - time0))
        analyze_suburb_parcel_data(metro_path, jobcenters_path, name, zoning_input)

    print('Finished')

def graph_results():

    # Get all zoning and lot data
    combined_jobcenter_zoning_data = pd.DataFrame()
    combined_ua_zoning_data = pd.DataFrame()
    combined_jobcenter_lotsize_data = pd.DataFrame()
    combined_ua_lotsize_data = pd.DataFrame()
    for name in ['Austin', 'Dallas', 'Houston']:

        def process_column_names(list_of_dfs):
            for df in list_of_dfs:
                df.columns = [str(col) + '-' + name for col in df.columns]

        # Get paths
        jc_lotsize_data = pd.read_csv(get_lotsize_path(name, 'jobcenter'), index_col = 0)
        ua_lotsize_data = pd.read_csv(get_lotsize_path(name, 'ua'), index_col = 0)
        jc_zoning_data = pd.read_csv(get_percent_zoned_path(name, 'jobcenter'), index_col = 0)
        ua_zoning_data = pd.read_csv(get_percent_zoned_path(name, 'ua'), index_col = 0)

        # Change column names
        process_column_names([jc_lotsize_data, ua_lotsize_data, jc_zoning_data, ua_zoning_data])

        # Combine
        combined_jobcenter_zoning_data = pd.concat([combined_jobcenter_zoning_data, jc_zoning_data], axis = 1)
        combined_ua_zoning_data = pd.concat([combined_ua_zoning_data, ua_zoning_data], axis = 1)
        combined_jobcenter_lotsize_data = pd.concat([combined_jobcenter_lotsize_data, jc_lotsize_data], axis = 1)
        combined_ua_lotsize_data = pd.concat([combined_ua_lotsize_data, ua_lotsize_data], axis = 1)

    for data, level in zip([combined_jobcenter_zoning_data, combined_ua_zoning_data], ['jobcenter', 'ua']):
        data['dist_from_center'] = data.index
        data = data.melt(var_name = 'zone_city', value_name = 'percent', id_vars=['dist_from_center'])
        data['City'] = data['zone_city'].apply(lambda x: str(x).split('-')[-1])
        data['Zone'] = data['zone_city'].apply(lambda x: str(x).split('-')[0])
        data = data.loc[data['dist_from_center'].apply(will_it_float)]
        data['dist_from_center'] = data['dist_from_center'].apply(lambda x: float(x))
        data = data.loc[(data['dist_from_center'] >= 10) & (data['dist_from_center'] <= 45)]
        data['percent'] = 100*data['percent']
        if level == 'ua':
            title = 'Land Use by Distance from City Center, All Suburbs'
        else:
            title = 'Land Use by Distance from City Center, Job Centers'
        plot = (ggplot(data, aes(x = 'dist_from_center', y = 'percent', fill = 'Zone'))
              + geom_col(position = 'dodge')
              + facet_wrap('~City')
              + labs(title = title, x = 'Distance from City Center', y = 'Percent of Land',
                     caption = 'Data based on Parcel Shapefiles from Counties throughout Texas'))
        path = 'Figures/Suburbs/all_{}_parcel_zone_rings.svg'.format(level)
        plot.save(path, width = 12, height = 8)

    for data, level in zip([combined_jobcenter_lotsize_data, combined_ua_lotsize_data], ['jobcenter', 'ua']):
        data['dist_from_center'] = data.index
        data = data.melt(var_name = 'zone_city', value_name = 'lotsize', id_vars=['dist_from_center'])
        data['City'] = data['zone_city'].apply(lambda x: str(x).split('-')[-1])
        data['Zone'] = data['zone_city'].apply(lambda x: str(x).split('-')[0])
        data = data.loc[data['dist_from_center'].apply(will_it_float)]
        data['dist_from_center'] = data['dist_from_center'].apply(lambda x: float(x))
        data = data.loc[(data['dist_from_center'] >= 10) & (data['dist_from_center'] <= 45)]
        if level == 'ua':
            title = 'Lotsize by Distance from City Center, All Suburbs'
        else:
            title = 'Lotsize by Distance from City Center, Job Centers'

        plot = (ggplot(data, aes(x = 'dist_from_center', y = 'lotsize', fill = 'City'))
              + geom_col(position = 'dodge')
              + facet_wrap('~Zone', scales = 'free')
              + labs(title = title, x = 'Distance from City Center', y = 'Lotsize (Square Feet)',
                     caption = 'Data based on Parcel Shapefiles from Counties throughout Texas'))
        path = 'Figures/Suburbs/all_{}_parcel_lotsize_rings.svg'.format(level)
        plot.save(path, width = 12, height = 8)


def calculate_area_and_distance_from_center(path, zoning_input, csv_path = None, write_flag = True):

    # Read file
    gdf = gpd.read_file(path)

    # Get area in square feet
    if 'sqft' not in gdf.columns:
        print('Calculating area')
        gdf = sf.get_area_in_units(gdf, scale=1, name='area_sqft')

    # Get distance from center of city, in miles currently
    if 'dist_to_cent' not in gdf.columns:
        print('Calculating distance')
        gdf['dist_to_cent'] = sf.calculate_dist_to_center(gdf, zoning_input)

    if write_flag:

        # Write data to csv
        if csv_path is not None:
            print('Writing csv')
            data = gdf[[col for col in gdf.columns if col != 'geometry']]
            data.to_csv(csv_path)

        # Write to shapefile
        gdf.to_file(path)

    return None

# At some point I ran:
    #for zoning_input, name in zip([austin_inputs, dallas_inputs, houston_inputs], ['austin', 'dallas', 'houston']):
        #for level in ['all', 'metro', 'jobcenter']:
            #path = get_cached_parcel_path(name, level)
            #csv_path = get_cached_parcel_path_csv(name, level)
            #calculate_area_and_distance_from_center(path, zoning_input, csv_path)
            #print('Finished with {} for {}, time is {}'.format(name, level, time.time() - time0))

def get_municipality_choropleth_path(name):
    return 'data/caches/suburbs/{}_municipality_calculations/municipality_calculations.shp'.format(name)

def calculate_municipalities(name, place_name, num_municipalities = 50, level = 'all', write_csv = True, to_simplify = None):
    """
    :param name: Name of the city, in all lowercase
    :param place_name: Name of the place (probably the name of the city with the first letter capitalized)
    :param num_municipalities: Number of municipalities to consider, defaults to 50 (will perform a nearest neighbor
    search and find the nearest ones).
    :param level: Level of the parcel data to use.
    :param write_csv: If true, will write the calculated parcel data, with municipalities, (minus the geometry) to a csv.
    :param to_simplify: Defaults to None. Should be a list of places to simplify the geometry for,
    because otherwise some intersections (i.e. with core municipalities) might take amount of time.
    :return: None, but caches the results.
    """

    # Get data
    data = pd.read_csv(get_cached_parcel_path_csv(name, level))
    data = data[['area_sqft', 'dist_to_cent']]
    parcel_data = gpd.read_file(get_cached_parcel_path(name, level))
    assert data.shape[0] == parcel_data.shape[0]
    parcel_data = parcel_data.merge(data, left_index = True, right_index = True)
    parcel_data['centroids'] = parcel_data['geometry'].centroid
    parcel_data = parcel_data.set_geometry('centroids')
    parcel_data['place'] = "Unincorporated"

    # Get place shapes
    texas_places = gpd.read_file(texas_places_path)
    places_spatial_index = texas_places.sindex
    city_polygon = texas_places.loc[texas_places['NAME'] == place_name, 'geometry']
    nearest_indexes = list(places_spatial_index.nearest(city_polygon.bounds.values[0], num_results = num_municipalities))
    nearest_places = texas_places.iloc[nearest_indexes]

    # Simplify nearest_places and get area
    if to_simplify is not None:
        nearest_places.loc[nearest_places['NAME'].isin(to_simplify), 'geometry'] = nearest_places.loc[nearest_places['NAME'].isin(to_simplify), 'geometry'].simplify(0.0001, preserve_topology = True)
    nearest_places = sf.get_area_in_units(nearest_places)

    # Now find which parcels are in which municipalities and calculate % sf, % mf, etc
    spatial_index = parcel_data.sindex
    counter = 0
    def process(row):
        """
        :param row: This is a row in the nearest_places dataframe!
        :return: The row but with eight more features (counts and total areas for sf, mf, other, and other res)
        """
        nonlocal parcel_data, counter
        print('Starting {}'.format(row["NAME"]))
        polygon = row['geometry']
        if row['area'] > 15: # more than 15 square miles --> fragment
            polygon_list = sf.fragment(polygon)
        else:
            polygon_list = [polygon]

        all_precise_intersections = set()
        for grid_piece in polygon_list:
            possible_intersections_index = list(spatial_index.intersection(grid_piece.bounds))
            possible_intersections = parcel_data.iloc[possible_intersections_index]
            precise_intersections = set(possible_intersections.loc[possible_intersections.intersects(grid_piece)].index.tolist())
            all_precise_intersections = all_precise_intersections.union(precise_intersections)
            print(len(all_precise_intersections))


        # Assign parcel data the precise intersections indexes
        parcel_data.loc[all_precise_intersections, 'place'] = row['NAME']

        # Calculate total area and number of points
        subset = parcel_data.loc[all_precise_intersections]
        total_area = subset.groupby(['broad_zone'])['area_sqft'].sum()/subset['area_sqft'].sum()
        total_area.index = [i + '-zone_pct' for i in total_area.index]
        num_points = subset.groupby(['broad_zone'])['area_sqft'].mean()
        num_points.index = [i + '-lotsize' for i in num_points.index]
        row = pd.concat([row, total_area, num_points])
        print('Finished with {} at {}'.format(row['NAME'], time.time() - time0, counter))

        counter += 1
        return row


    nearest_places = nearest_places.apply(process, axis = 1)

    outfile_path = get_municipality_choropleth_path(name)
    nearest_places.to_file(outfile_path)

    print(parcel_data.loc[parcel_data['place'] != 'Unincorporated'].shape[0]/parcel_data.shape[0])

    if write_csv:

        data = parcel_data[[col for col in parcel_data.columns if col != 'geometry']]
        path = get_cached_parcel_path_csv(name, level)
        print('Writing csv')
        data.to_csv(path)
        print('Finished, time is {}'.format(time.time() - time0))

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

    # Read data
    def land_use_by_municipality_path(name):
        return 'data/caches/suburbs/{}_land_use_by_municipality.csv'.format(name)

    path = get_cached_parcel_path_csv(name, 'all')
    data = pd.read_csv(path)

    # Drop duplicates by centroids - note that these centroids are strings and have not been parsed as points yet
    data = data.drop_duplicates(subset = 'centroids', keep = 'first')

    # Calculate percent of area zoned
    zone_areas = data.groupby(['broad_zone', 'place'])['area_sqft'].sum()
    municipality_areas = data.groupby(['place'])['area_sqft'].sum()
    final_data = zone_areas.divide(municipality_areas)
    final_data = final_data.reset_index()
    final_data.columns = ['broad_zone', 'place', 'Percent of Land Used']
    print(final_data)
    final_data.to_csv(land_use_by_municipality_path(name))



if __name__ == '__main__':

    #analyze_transportation_networks(['austin', 'dallas', 'houston'], [austin_inputs, dallas_inputs, houston_inputs])
    #plot_municipality_choropleth('austin', austin_inputs, austin_job_centers)
    #plot_municipality_choropleth('dallas', dallas_inputs, dallas_job_centers, to_exclude = ['Combine'])
    #plot_municipality_choropleth('houston', houston_inputs, houston_job_centers)

    analyze_land_use_by_metro('austin')
    analyze_land_use_by_metro('dallas')
    analyze_land_use_by_metro('houston')

