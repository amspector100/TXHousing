"""Graphs which mostly rely on parcel data"""
import time
import numpy as np
import pandas as pd
import geopandas as gpd
from .. import utilities
from ..data_processing import zoning, parcel
from plotnine import *

def plot_singlefamily_lotsizes(calculate = False, cache_path = 'data/caches/Cached_Lotsize_Data_v2.csv',
                               save_path = 'Figures/Bucket 2/sf_lotsizes.svg', width = 10, height = 8):
    """Plots the actual average single family lotsize by distance from city center in the TX Triangle, calculated from
    parcel data.

    :param calculate: If true, recalculate the averages. If not, will try to read the data from the cache_path.
    :param cache_path: The location at which to save a csv of cached averages.
    :param save_path: The location at which to save the final graph.
    :param width: The width of the final graph, in inches or pixels (depends on file format)
    :param height: The height of the final graph, in inches or pixels (depends on the file format)"""

    time0 = time.time()

    if calculate:

        maximum = 10

        # Initialize result
        def get_rings_of_parcels(parcel_data, feature, zone_dictionary, zoning_input, cityname):

            def parse_broad_zone(text):
                for key in zone_dictionary:
                    if text in zone_dictionary[key]:
                        return key
                return 'Other'

            parcel_data['broad_zone'] = parcel_data[feature].apply(parse_broad_zone)
            parcel_data = parcel_data.loc[parcel_data['broad_zone'].isin(['Single Family', 'Multifamily'])]

            # Get the area in square feet and then get the mean lot size by distance from city center
            parcel_data = utilities.measurements.get_area_in_units(parcel_data, scale=1, name='area')

            # For efficiency, create centroids and use these
            parcel_data['centroids'] = parcel_data['geometry'].centroid
            parcel_data = parcel_data.set_geometry('centroids')

            # Initialize result
            result = pd.DataFrame()

            for zone in ['Single Family', 'Multifamily']:
                cityzone = str(cityname) + '-' + zone
                print(cityzone)
                result[cityzone] = utilities.measurements.points_intersect_rings(parcel_data.loc[parcel_data['broad_zone'] == zone],
                                                                                 lat = zoning_input.lat,
                                                                                 long = zoning_input.long,
                                                                                 factor = 'area',
                                                                                 step = 1,
                                                                                 categorical = False,
                                                                                 by = 'mean',
                                                                                 per_square_mile = False,
                                                                                 geometry_column = 'centroids',
                                                                                 maximum = maximum)

            return result

        # Work on Houston data ---------------------------------------------------
        print('Starting to work on Houston data, time is {}'.format(time.time() - time0))
        houston_parcels = parcel.process_houston_parcel_data()
        print('Finished reading Houston parcel data, time is {}'.format(time.time() - time0))
        houston_parcels = houston_parcels.replace([np.inf, -np.inf], np.nan)
        houston_parcels = houston_parcels.dropna(subset = ['USE_CODE'], how = 'all')
        houston_parcels['USE_CODE'] = houston_parcels['USE_CODE'].astype(str)
        houston_feature = 'USE_CODE'

        # Here, 1006 refers to condominiums, 1007 refers to townhomes
        houston_dictionary = parcel.state_sptbcode_dictionary
        houston_result = get_rings_of_parcels(houston_parcels, houston_feature, houston_dictionary,
                                              zoning.houston_inputs, 'Houston')

        # Work on Austin data ---------------------------------------------------
        print('Starting to work on Austin data, time is {}'.format(time.time() - time0))
        austin_parcels = parcel.process_austin_parcel_data()
        print('Finished reading austin parcel data, time is {}'.format(time.time() - time0))
        austin_feature = 'basezone'
        # Note SF-4B refers to condominium, and SF-6 refers to 'townhouse and condominium'
        austin_dictionary = {'Single Family':['SF-1', 'SF-2', 'SF-3', 'SF-4A', 'SF-5'],
                             'Multifamily':['SF-4B', 'SF-6', 'MF-1', 'MF-2', 'MF-3', 'MF-4', 'MF-5', 'MF-6']}

        austin_result = get_rings_of_parcels(austin_parcels, austin_feature, austin_dictionary,
                                             zoning.austin_inputs, 'Austin')

        # Work on Dallas data ---------------------------------------------------
        print('Starting to work on Dallas data, time is {}'.format(time.time() - time0))
        dallas_parcels = gpd.read_file(parcel.dallas_parcel_data_path_2016)
        print('Finished reading dallas parcel data, time is {}'.format(time.time() - time0))
        dallas_feature = 'sptbcode'
        # Note A12 are townhouses, A13 are condominiums
        dallas_dictionary = parcel.dallas_sptb_dictionary
        dallas_result = get_rings_of_parcels(dallas_parcels, dallas_feature, dallas_dictionary, zoning.dallas_inputs,
                                             'Dallas')

        all_results = pd.concat([austin_result, dallas_result, houston_result], axis = 1)
        all_results.to_csv(cache_path)

    all_results = pd.read_csv(cache_path)
    all_results = all_results.melt(var_name = 'Zone', value_name = 'avg_lot_size', id_vars = ['dist_to_center'])
    all_results['City'] = all_results['Zone'].apply(lambda x: x.split('-')[0])
    all_results['Zone'] = all_results['Zone'].apply(lambda x: x.split('-')[1])

    # Subset and graph
    all_results.index = all_results['dist_to_center']
    sf_results = utilities.measurements.order_radii(all_results.loc[all_results['Zone'] == 'Single Family'])
    sf_results['dist_to_center'] = sf_results.index
    sflotplot = (ggplot(sf_results, aes(x = 'dist_to_center', y = 'avg_lot_size', fill = 'City'))
                    + geom_col(position = 'dodge', width = 0.7)
                    + labs(x = 'Distance from Center of City (Miles', y = 'Average Lot Size (Square Feet)',
                           title = 'Average Lot Sizes by Distance from City Center, in Austin, Dallas, and Houston',
                           caption = 'Based on Parcel Data provided by Austin, Dallas, and Harris County.')
                    + theme_bw())
    sflotplot.save(save_path, width = width, height = height)

def get_parcel_feature_outfile(name):
    """Path for caches used in plot_percent_undeveloped"""
    return 'data/caches/{}_municipal_parcel_features.csv'.format(name)


def plot_percent_undeveloped():

    from suburbs import process_lot_descriptions

    warning_flag = True

    if calculate and warning_flag:
        raise Warning(
            'You just tried to recalculate the percent of land that is undeveloped in each city. The problem is that'
            'redoing this will actually add extra columns on and mess up the CSVs. On this run, I threw an exception.'
            'If you are totally sure you want to do this, change "warning_flag" to False in line 214 in misc_calcs.py.'
            'But make sure it is really what you want.')
    elif calculate:

        def fill_negatives_with_zeroes(a_series):
            a_series[a_series < 0] = 0
            return a_series

        # The drill is: read in data, process broad_zone, then get percent undeveloped --

        # Austin
        austin_data = pd.read_csv(get_parcel_feature_outfile(
            'austin'))  # [['val_per_sqft', 'dist_to_center', 'appraised', 'area', 'far', 'lu_desc']]
        austin_dictionary = {'Single Family': ['Single Family', 'Duplexes', 'Large-lot Single Family'],
                             'Multifamily': ['Apartment/Condo', 'Three/Fourplex', 'Group Quarters', 'Mixed Use']}
        austin_data['lu_desc'] = austin_data['lu_desc'].astype(str)
        austin_data['broad_zone'] = process_lot_descriptions(austin_data, 'lu_desc', austin_dictionary)
        austin_data['percent_undeveloped'] = fill_negatives_with_zeroes(1 - austin_data['far'])
        austin_data.drop('Unnamed: 0', axis=1, inplace=True)
        austin_data.to_csv(get_parcel_feature_outfile('austin'))

        # Dallas
        dallas_data = pd.read_csv(get_parcel_feature_outfile(
            'dallas'))  # [['Acct', 'val_per_sqft', 'dist_to_center', 'TOT_VAL', 'area', 'SPTD_CODE']]
        dallas_data['broad_zone'] = process_lot_descriptions(dallas_data, 'SPTD_CODE', dallas_sptb_dictionary)
        dallas_extra_data = pd.read_csv(dallas_county_res_path)
        dallas_extra_data = dallas_extra_data.drop_duplicates(subset='ACCOUNT_NUM', keep='first')
        dallas_data = dallas_data.merge(dallas_extra_data, left_on='Acct', right_on='ACCOUNT_NUM', how='left')
        dallas_data['percent_undeveloped'] = fill_negatives_with_zeroes(
            1 - dallas_data['TOT_MAIN_SF'].divide(dallas_data['area']))
        dallas_data.loc[dallas_data['broad_zone'] == "Other", 'percent_undeveloped'] = float("NaN")
        dallas_data.drop('Unnamed: 0', axis=1, inplace=True)
        dallas_data.to_csv(get_parcel_feature_outfile('dallas'))

        # Houston
        houston_data = pd.read_csv(get_parcel_feature_outfile(
            'houston'))  # [['val_per_sqft', 'dist_to_center', 'TOTAL_APPRAISED_VALUE', 'area', 'STATE_CLASS']]
        houston_data['broad_zone'] = process_lot_descriptions(houston_data, 'STATE_CLASS',
                                                              state_sptbcode_dictionary)
        houston_extra_data = pd.read_csv(harris_parcel_building_res_path_2018, sep='\t', header=None,
                                         encoding='Latin-1')
        houston_extra_data.columns = houston_building_res_columns
        houston_extra_data['ACCOUNT'] = houston_extra_data['ACCOUNT'].astype(float)
        houston_extra_data = houston_extra_data.drop_duplicates(subset='ACCOUNT', keep='first')
        houston_data = houston_data.merge(houston_extra_data, left_on='HCAD_NUM', right_on='ACCOUNT')
        houston_data['percent_undeveloped'] = fill_negatives_with_zeroes(
            1 - houston_data['BASE_AREA'].divide(houston_data['area']))
        houston_data.loc[houston_data['broad_zone'] == "Other", 'percent_undeveloped'] = float("NaN")
        houston_data.drop('Unnamed: 0', axis=1, inplace=True)
        houston_data.to_csv(get_parcel_feature_outfile('houston'))

    # Combine data, cache result, and graph
    austin_data = pd.read_csv(get_parcel_feature_outfile('austin'))[
        ['dist_to_center', 'broad_zone', 'percent_undeveloped', 'far', 'area']]
    austin_data['developed_sqft'] = pd.concat(
        [austin_data['area'].multiply(austin_data['far']), austin_data['area']], axis=1).min(axis=1)
    austin_data['City'] = 'Austin'
    dallas_data = pd.read_csv(get_parcel_feature_outfile('dallas'))[
        ['dist_to_center', 'broad_zone', 'percent_undeveloped', 'TOT_MAIN_SF', 'area']]
    dallas_data['City'] = 'Dallas'
    dallas_data['developed_sqft'] = pd.concat([dallas_data['TOT_MAIN_SF'], dallas_data['area']], axis=1).min(axis=1)
    houston_data = pd.read_csv(get_parcel_feature_outfile('houston'))[
        ['dist_to_center', 'broad_zone', 'percent_undeveloped', 'BASE_AREA', 'area']]
    houston_data['developed_sqft'] = pd.concat([houston_data['BASE_AREA'], houston_data['area']], axis=1).min(
        axis=1)
    houston_data['City'] = 'Houston'

    # Combined
    all_data = pd.concat([austin_data, dallas_data, houston_data])
    all_data['smoothed_dist_to_center'] = all_data['dist_to_center'].apply(np.ceil)

    # Medians
    medians = (100 * (all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])[
        'percent_undeveloped'].median())).reset_index()
    medians = medians.rename(columns={'percent_undeveloped': 'percent_undeveloped_median'})

    # Simple means
    simple_means = (100 * (all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])[
        'percent_undeveloped'].mean())).reset_index()
    simple_means = simple_means.rename(columns={'percent_undeveloped': 'percent_undeveloped_simple_mean'})

    # Weighted means
    total_areas = all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])['area'].sum()
    total_developed_area = (100 * (
    all_data.groupby(['smoothed_dist_to_center', 'broad_zone', 'City'])['developed_sqft'].sum(skipna=True)))
    weighted_means = total_developed_area.divide(total_areas)

    def replace_high_values(x):
        if x > 100:
            return 100
        else:
            return x

    weighted_means = weighted_means.apply(replace_high_values).reset_index().rename(
        columns={0: 'percent_undeveloped_weighted_average'})

    all_results = reduce(
        lambda left, right: pd.merge(left, right, on=['broad_zone', 'smoothed_dist_to_center', 'City'], how='outer',
                                     sort=False),
        [weighted_means, simple_means, medians])
    all_results.to_csv('shared_data/percent_undeveloped.csv')

    # Now graph
    all_results = all_results.loc[all_results['broad_zone'] == 'Single Family']
    all_results = all_results.loc[all_results['smoothed_dist_to_center'] < 11]
    p = (ggplot(all_results,
                aes(x='smoothed_dist_to_center', y='percent_undeveloped_simple_mean', fill='City', group='City'))
         + geom_col(position='dodge')
         + labs(title='Undeveloped Area of Single Family Lots in Texas Triangle',
                x='Distance from City Center',
                y='Simple Mean of Undeveloped Percentage of SF Lots'))
    print(p)
    p.save('Figures/Bucket 2/percent_undeveloped.svg', width=10, height=8)