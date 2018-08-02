from inputs import *
from helpers import *
import spatial_functions as sf
from tqdm import tqdm
from functools import reduce

import time
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import folium
from plotnine import *


special_setbacks_and_lots = False
calc_percent_residential_in_block_groups = False
calc_landmarks_in_houston = False
calc_parking_cost = False
calc_percent_used = True
validate_parcel_results = False
calc_overlays = False

calculate = False

time0 = time.time()

# Special setback and minimum lot zone calculations ----------------------------- (Houston) ----------------------------------------------------- Special setback and minimum lot zone calculations
if special_setbacks_and_lots:

    print('Running some numbers on special setbacks and lots')
    special_setbacks = gpd.read_file(houston_spec_min_setbacks)
    special_lots = gpd.read_file(houston_spec_min_lots)

    # Find areas
    special_setbacks = sf.get_area_in_units(special_setbacks) # Remember, this defaults to miles
    special_lots = sf.get_area_in_units(special_lots)
    print('Special setbacks area', special_setbacks['area'].sum())
    print('Special lots area', special_lots['area'].sum())

    # Find population in the particular areas - B01001e1 is the column which counts population in he block data
    block_data = sf.get_block_geodata(['X01_AGE_AND_SEX'], cities = 'Houston')
    special_lots = sf.get_all_averages_by_area(block_data, special_lots, fillna = 0) # Defaults to getting pop data
    special_setbacks = sf.get_all_averages_by_area(block_data, special_setbacks, fillna = 0)

    print('Special lots population', special_lots['B01001e1'].sum())
    print('Special setbacks population', special_setbacks['B01001e1'].sum())


# Calculate percent of land that is residential for block groups -------------------- (Austin/Dallas) ---------------------------------------------------- Calculate percent of land that is residential

if calc_percent_residential_in_block_groups:
    # Get block and zoning data
    block_data = sf.get_block_geodata(['X19_INCOME'], cities = ['Austin', 'Dallas'])
    austin_block_data = block_data['Austin']
    dallas_block_data = block_data['Dallas']

    # Get zone data and simplify for performance. Note zone data is a geoseries not a geodataframe.
    dallas_zones = process_zoning_shapefile(dallas_inputs, broaden = True).to_crs({'init':'epsg:4326'})
    dallas_zones = dallas_zones.loc[dallas_zones['broad_zone'].isin(['Single Family', 'Multifamily'])]
    dallas_zones = dallas_zones['geometry']
    austin_zones = sf.process_geometry(process_zoning_shapefile(austin_inputs, broaden = True))
    austin_zones = austin_zones.loc[austin_zones['broad_zone'].isin(['Single Family', 'Multifamily'])]
    austin_zones = austin_zones['geometry']


    # Now find intersections
    for cityname, block_data, zone_data in zip(['Dallas', 'Austin'], [dallas_block_data, austin_block_data], [dallas_zones, austin_zones]):
        print('Starting to work on {}'.format(cityname))
        spatial_index = zone_data.sindex
        block_data['area'] = block_data['geometry'].area

        def calc_percent_res(polygon):
            poly_area = polygon.area
            possible_intersections_index = list(spatial_index.intersection(polygon.bounds))
            possible_intersections = zone_data.iloc[possible_intersections_index]
            precise_intersections = possible_intersections.intersection(polygon)
            res_area = precise_intersections.area.sum()
            percent_residential = res_area/poly_area
            return percent_residential

        block_data['percent_residential'] = block_data['geometry'].apply(calc_percent_res)

    results = pd.concat([austin_block_data['percent_residential'], dallas_block_data['percent_residential']])
    results.to_csv('shared_data/bg_percent_residential.csv')

# Calculate population of nat hist districts in Houston as well as
if calc_landmarks_in_houston:

    # Get points and polygons
    houston_landmarks = gpd.read_file(houston_historic_landmarks_path).to_crs({'init':'epsg:4326'})
    tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
    tx_hd_data = gpd.read_file(tx_hd_path)
    houston_nat_districts = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston'].to_crs({'init':'epsg:4326'})
    houston_block_data = sf.get_block_geodata(['X01_AGE_AND_SEX'], cities = 'Houston')

    # Calculate number of points in polygons
    num_points = sf.points_over_area(houston_landmarks, houston_nat_districts, normalize_by_area = False)
    population = sf.get_all_averages_by_area(houston_block_data, houston_nat_districts, fillna=0,
                                             account_method='water')['B01001e1'].sum()

    print('Number of Historic Landmarks in NRHDS is {}, Population of NRHDs is {}'.format(num_points, population))

# Calculate parking cost
if calc_parking_cost:


    if calculate:

        def subset_to_core(name, path, zoning_input, value_feature, area_feature = None, merge_path = None, left_on = None, right_on = None, parcel_data = None, write_file = True):

            if parcel_data is None:
                # Read parcel data
                print('Reading, time is {}'.format(time.time() - time0))
                parcel_data = gpd.read_file(path)

                # Merge
                if merge_path is not None:
                    if left_on is None or right_on is None:
                        raise ValueError('When merge_path is not None, left_on and right_on must also not be None in subset_to_core function')
                    extra_data = pd.read_csv(merge_path)
                    parcel_data[left_on] = parcel_data[left_on].astype(str)
                    extra_data[right_on] = extra_data[right_on].astype(str)
                    extra_data = extra_data.drop_duplicates(subset = right_on, keep = 'first')
                    parcel_data = parcel_data.merge(extra_data, left_on = left_on, right_on = right_on, how = 'left')

            # Transform and drop duplicates
            parcel_data['centroids_string'] = parcel_data['geometry'].centroid.astype(str)
            parcel_data = parcel_data.drop_duplicates(subset = 'centroids_string', keep = 'first')

            if parcel_data.crs != {'init':'epsg:4326'} and parcel_data.crs is not None:
                parcel_data = parcel_data.to_crs({'init':'epsg:4326'})
            elif parcel_data.crs is None:
                print('No CRS for parcel data, assuming it is in lat long coordinates.')

            # Get dist to center
            parcel_data['dist_to_center'] = sf.calculate_dist_to_center(parcel_data, zoning_input, drop_centroids=False)

            # Subset
            core_parcel_data = parcel_data.loc[(parcel_data[value_feature] != 0) & (parcel_data[value_feature].notnull())]

            # Get area
            if area_feature is None:
                print('Getting area, time is {}'.format(time.time() - time0))
                core_parcel_data = sf.get_area_in_units(core_parcel_data, name = 'area', scale = 1)
            else:
                core_parcel_data = core_parcel_data.rename(columns = {area_feature: 'area'})

            # Get val per sqft
            core_parcel_data['val_per_sqft'] = core_parcel_data[value_feature].divide(core_parcel_data['area'])

            # Write to csv (cache)
            if write_file:
                core_parcel_data[[col for col in core_parcel_data.columns if col not in ['centroids', 'geometry']]].to_csv(get_parcel_feature_outfile(name))
            else:
                print('write_file is False, so not caching data')

            return core_parcel_data

        # Austin
        core_austin_data = subset_to_core('austin', austin_parcel_path, austin_inputs, 'appraised') # Best data

        # Checked, Dallas's 'area_feet' column matches with my own area calculations
        core_dallas_data = subset_to_core('dallas', dallas_county_parcel_path, dallas_inputs, value_feature = 'TOT_VAL', # Mediocre data
                                          merge_path = dallas_county_appraisal_path, left_on = 'Acct', right_on = 'ACCOUNT_NUM')

        # Houston
        houston_parcel_data = process_houston_parcel_data(feature_files=[harris_parcel_land_path_2018, harris_parcel_appraisal_path_2018], # Mediocre data
                                    feature_columns_list=[houston_land_columns, houston_appraisal_columns], county_level=False)
        core_houston_data = subset_to_core('houston', None, houston_inputs, 'TOTAL_APPRAISED_VALUE', parcel_data = houston_parcel_data)

    # Core radius
    core_radius = 1

    # Now read csvs
    austin_data = pd.read_csv(get_parcel_feature_outfile('austin'))[['val_per_sqft', 'dist_to_center', 'appraised', 'area']]
    austin_data = austin_data.astype(float).rename(columns = {'appraised':'value'})
    austin_data['City'] = 'Austin'
    dallas_data = pd.read_csv(get_parcel_feature_outfile('dallas'))[['val_per_sqft', 'dist_to_center', 'TOT_VAL', 'area']]
    dallas_data = dallas_data.astype(float).rename(columns = {'TOT_VAL':'value'})
    dallas_data['City'] = 'Dallas'
    houston_data = pd.read_csv(get_parcel_feature_outfile('houston'))[['val_per_sqft', 'dist_to_center', 'TOTAL_APPRAISED_VALUE', 'area']]
    houston_data = houston_data.astype(float).rename(columns = {'TOTAL_APPRAISED_VALUE':'value'})
    houston_data['City'] = 'Houston'

    all_data = pd.concat([austin_data, dallas_data, houston_data])
    all_data = all_data.loc[all_data['dist_to_center'] <= core_radius]
    print('Median:')
    print(all_data.groupby(['City'])['val_per_sqft'].median())
    print('Simple mean:')
    print(all_data.groupby(['City'])['val_per_sqft'].mean())
    print('Weighted mean:')
    print(all_data.groupby(['City'])['value'].sum().divide(all_data.groupby(['City'])['area'].sum()))

    # Get rid of the outliers
    maximum = all_data['val_per_sqft'].quantile(.95)
    all_data.loc[all_data['val_per_sqft'] > maximum, 'val_per_sqft'] = maximum

    p = (ggplot(all_data, aes(x = 'val_per_sqft', fill = 'City', group = 'City'))
          + geom_histogram()
          + labs(title = 'Property Values within 1 Mile of City Center in Texas Triangle',
                 x = 'Value per Square Foot',
                 y = 'Number of Parcels')
          + facet_wrap('~City'))
    print(p)
    p.save('Figures/property_value_histogram.svg', width = 10, height = 8)

if calc_percent_used:

    from suburbs import process_lot_descriptions

    warning_flag = True

    if calculate and warning_flag:
        raise Warning('You just tried to recalculate the percent of land that is undeveloped in each city. The problem is that'
                      'redoing this will actually add extra columns on and mess up the CSVs. On this run, I threw an exception.'
                      'If you are totally sure you want to do this, change "warning_flag" to False in line 214 in misc_calcs.py.'
                      'But make sure it is really what you want.')
    elif calculate:

        def fill_negatives_with_zeroes(a_series):
            a_series[a_series < 0] = 0
            return a_series

        # The drill is: read in data, process broad_zone, then get percent undeveloped --

        # Austin
        austin_data = pd.read_csv(get_parcel_feature_outfile('austin'))#[['val_per_sqft', 'dist_to_center', 'appraised', 'area', 'far', 'lu_desc']]
        austin_dictionary = {'Single Family':['Single Family', 'Duplexes', 'Large-lot Single Family'], 'Multifamily':['Apartment/Condo', 'Three/Fourplex', 'Group Quarters', 'Mixed Use']}
        austin_data['lu_desc'] = austin_data['lu_desc'].astype(str)
        austin_data['broad_zone'] = process_lot_descriptions(austin_data, 'lu_desc', austin_dictionary)
        austin_data['percent_undeveloped'] = fill_negatives_with_zeroes(1 - austin_data['far'])
        austin_data.drop('Unnamed: 0', axis = 1, inplace = True)
        austin_data.to_csv(get_parcel_feature_outfile('austin'))

        # Dallas
        dallas_data = pd.read_csv(get_parcel_feature_outfile('dallas'))#[['Acct', 'val_per_sqft', 'dist_to_center', 'TOT_VAL', 'area', 'SPTD_CODE']]
        dallas_data['broad_zone'] = process_lot_descriptions(dallas_data, 'SPTD_CODE', dallas_sptb_dictionary)
        dallas_extra_data = pd.read_csv(dallas_county_res_path)
        dallas_extra_data = dallas_extra_data.drop_duplicates(subset = 'ACCOUNT_NUM', keep = 'first')
        dallas_data = dallas_data.merge(dallas_extra_data, left_on = 'Acct', right_on = 'ACCOUNT_NUM', how = 'left')
        dallas_data['percent_undeveloped'] = fill_negatives_with_zeroes(1 - dallas_data['TOT_MAIN_SF'].divide(dallas_data['area']))
        dallas_data.loc[dallas_data['broad_zone'] == "Other", 'percent_undeveloped'] = float("NaN")
        dallas_data.drop('Unnamed: 0', axis = 1, inplace = True)
        dallas_data.to_csv(get_parcel_feature_outfile('dallas'))

        # Houston
        houston_data = pd.read_csv(get_parcel_feature_outfile('houston'))#[['val_per_sqft', 'dist_to_center', 'TOTAL_APPRAISED_VALUE', 'area', 'STATE_CLASS']]
        houston_data['broad_zone'] = process_lot_descriptions(houston_data, 'STATE_CLASS', state_sptbcode_dictionary)
        houston_extra_data = pd.read_csv(harris_parcel_building_res_path_2018, sep='\t', header=None, encoding='Latin-1')
        houston_extra_data.columns = houston_building_res_columns
        houston_extra_data['ACCOUNT'] = houston_extra_data['ACCOUNT'].astype(float)
        houston_extra_data = houston_extra_data.drop_duplicates(subset = 'ACCOUNT', keep = 'first')
        houston_data = houston_data.merge(houston_extra_data, left_on = 'HCAD_NUM', right_on = 'ACCOUNT')
        houston_data['percent_undeveloped'] = fill_negatives_with_zeroes(1 - houston_data['BASE_AREA'].divide(houston_data['area']))
        houston_data.loc[houston_data['broad_zone'] == "Other", 'percent_undeveloped'] = float("NaN")
        houston_data.drop('Unnamed: 0', axis = 1, inplace = True)
        houston_data.to_csv(get_parcel_feature_outfile('houston'))

    # Combine data, cache result, and graph
    austin_data = pd.read_csv(get_parcel_feature_outfile('austin'))[['dist_to_center', 'broad_zone', 'percent_undeveloped', 'far', 'area']]
    austin_data['developed_sqft'] = pd.concat([austin_data['area'].multiply(austin_data['far']), austin_data['area']], axis = 1).min(axis = 1)
    austin_data['City'] = 'Austin'
    dallas_data = pd.read_csv(get_parcel_feature_outfile('dallas'))[['dist_to_center', 'broad_zone', 'percent_undeveloped', 'TOT_MAIN_SF', 'area']]
    dallas_data['City'] = 'Dallas'
    dallas_data['developed_sqft'] = pd.concat([dallas_data['TOT_MAIN_SF'], dallas_data['area']], axis = 1).min(axis = 1)
    houston_data = pd.read_csv(get_parcel_feature_outfile('houston'))[['dist_to_center', 'broad_zone', 'percent_undeveloped', 'BASE_AREA', 'area']]
    houston_data['developed_sqft'] = pd.concat([houston_data['BASE_AREA'], houston_data['area']], axis = 1).min(axis = 1)
    houston_data['City'] = 'Houston'

    # Combined
    all_data = pd.concat([austin_data, dallas_data, houston_data])
    all_data['smoothed_dist_to_center'] = all_data['dist_to_center'].apply(np.ceil)

    # Medians
    medians = (100*(all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])['percent_undeveloped'].median())).reset_index()
    medians = medians.rename(columns = {'percent_undeveloped':'percent_undeveloped_median'})

    # Simple means
    simple_means = (100*(all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])['percent_undeveloped'].mean())).reset_index()
    simple_means = simple_means.rename(columns = {'percent_undeveloped':'percent_undeveloped_simple_mean'})

    # Weighted means
    total_areas = all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])['area'].sum()
    total_developed_area = (100*(all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])['developed_sqft'].sum(skipna = True)))
    weighted_means = total_developed_area.divide(total_areas)

    def replace_high_values(x):
        if x > 100:
            return 100
        else:
            return x

    weighted_means = weighted_means.apply(replace_high_values).reset_index().rename(columns = {0:'percent_undeveloped_weighted_average'})

    all_results = reduce(lambda left, right: pd.merge(left, right, on = ['broad_zone', 'smoothed_dist_to_center', 'City'], how = 'outer', sort = False),
                         [weighted_means, simple_means, medians])
    all_results.to_csv('shared_data/percent_undeveloped.csv')

    # Now graph
    all_results = all_results.loc[all_results['broad_zone'] == 'Single Family']
    all_results = all_results.loc[all_results['smoothed_dist_to_center'] < 11]
    p = (ggplot(all_results, aes(x = 'smoothed_dist_to_center', y = 'percent_undeveloped_simple_mean', fill = 'City', group = 'City'))
         + geom_col(position = 'dodge')
         + labs(title = 'Undeveloped Area of Single Family Lots in Texas Triangle',
                x = 'Distance from City Center',
                y = 'Simple Mean of Undeveloped Percentage of SF Lots'))
    print(p)
    p.save('Figures/Bucket 2/percent_undeveloped.svg', width = 10, height = 8)

if validate_parcel_results:

    from suburbs import all_dallas_zoning_path_csv

    zdata = pd.read_csv('data/caches/suburbs/dallas_zoning_use_by_municipality.csv', index_col = 0)
    pdata = pd.read_csv('data/caches/suburbs/{}_land_use_by_municipality.csv'.format('dallas'), index_col = 0)
    alldata = pdata.merge(zdata, on = ['broad_zone', 'place'], how = 'inner')
    alldata = alldata.loc[alldata['broad_zone'] == 'Single Family']
    difference = alldata['Percent of Land Used'] - alldata['Percent of Land Zoned As']
    print(difference.mean())
    print(difference.apply(abs).mean())

if calc_overlays:

    austin_zones = process_zoning_shapefile(austin_inputs)
    def get_overlay(row):
        row['overlay'] = str(row['zoning_zty']).replace(str(row['base_zone']), '')
        return row
    austin_zones = austin_zones.apply(get_overlay, axis = 1)
    austin_zones.crs = {'init':'epsg:4326'}
    austin_zones = sf.get_area_in_units(austin_zones)
    austin_zones['dist_to_center'] = sf.calculate_dist_to_center(austin_zones, austin_inputs)
    austin_zones['dist_to_center'] = austin_zones['dist_to_center'].apply(lambda x: 2*np.ceil(x/2))
    austin_zones['num_overlays'] = austin_zones['overlay'].apply(lambda x: str(x).count('-'))
    result = austin_zones.groupby(['dist_to_center', 'num_overlays'])['area'].count()
    divisor = austin_zones.groupby(['dist_to_center'])['area'].count()
    print(result, divisor)
    result = result.divide(divisor)
    print(result)