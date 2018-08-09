"""Graphs which mostly rely on municipal parcel data. Note that all of these graphs rely on the cached parcel data in
 csv format, not the actual parcel data, which is processed and cached in the data_processing package."""


import time
import numpy as np
import pandas as pd
import geopandas as gpd
from .. import utilities
from ..data_processing import zoning, parcel
from plotnine import *
from functools import reduce

def plot_singlefamily_lotsizes(save_path = 'Figures/Zoning/sf_lotsizes.svg', width = 10, height = 8):
    """ Plots average lotsizes of single family homes in Houston conditional on distance from city center. This uses
    cached municipal parcel data. """

    # Read data
    austin_path = parcel.get_cached_municipal_parcel_path('austin')
    austin_parcels = pd.read_csv(austin_path)[['dist_to_center', 'area_sqft', 'broad_zone', 'place']]
    austin_parcels = austin_parcels.loc[(austin_parcels['place'] == 'Austin')
                                        & (austin_parcels['broad_zone'] == 'Single Family')]

    dallas_path = parcel.get_cached_municipal_parcel_path('dallas')
    dallas_parcels = pd.read_csv(dallas_path)[['dist_to_center', 'area_sqft', 'broad_zone', 'place']]
    dallas_parcels = dallas_parcels.loc[(dallas_parcels['place'] == 'Dallas')
                                        & (dallas_parcels['broad_zone'] == 'Single Family')]

    houston_path = parcel.get_cached_municipal_parcel_path('houston')
    houston_parcels = pd.read_csv(houston_path)[['dist_to_center', 'area_sqft', 'broad_zone', 'place']]
    houston_parcels = houston_parcels.loc[(houston_parcels['place'] == 'Houston')
                                          & (houston_parcels['broad_zone'] == 'Single Family')]

    # Combine, smooth dist_to_center
    all_parcels = pd.concat([austin_parcels, dallas_parcels, houston_parcels], axis = 0, ignore_index = True)
    all_parcels['dist_to_center'] = all_parcels['dist_to_center'].astype(float).apply(np.ceil)
    all_parcels['dist_to_center'] = all_parcels['dist_to_center'].apply(lambda x: x if x <= 10 else '10+')

    # Calculate final graphed result
    result = all_parcels.groupby(['place', 'dist_to_center', 'broad_zone'])['area_sqft'].mean()
    result = result.unstack().reset_index()
    print(result.columns)
    sflotplot = (ggplot(result, aes(x = 'dist_to_center', y = 'Single Family', fill = 'place'))
                    + geom_col(position = 'dodge', width = 0.7)
                    + labs(x = 'Distance from Center of City (Miles', y = 'Average Lot Size (Square Feet)',
                           title = 'Average Lot Sizes by Distance from City Center, in Austin, Dallas, and Houston',
                           caption = 'Based on Parcel Data provided by Austin, Dallas, and Harris County.')
                    + theme_bw())
    sflotplot.save(save_path, width = width, height = height)



def plot_percent_undeveloped(save_path = 'Figures/Zoning/percent_undeveloped.svg'):
    """ Calculates the percent of land which is undeveloped (using base_area features and area calculations) conditional
    on distance from city center and broad_zone. Based on cached municipal parcel data."""


    # Some buildings/parcels span multiple zones, making their % undeveloped appear to be negative. We just substitute
    # these negative values with a value of 0, because most of these buildings are in very dense areas with extremely
    # low setback requirements anyway. This is a limitation on the accuracy of the results, of course. 
    def fill_negatives_with_zeroes(a_series):
        a_series[a_series < 0] = 0
        return a_series

    # Read data and calculate percent undeveloped as well as the number of developed square feet
    austin_path = parcel.get_cached_municipal_parcel_path('austin')
    austin_parcels = pd.read_csv(austin_path)[['dist_to_center', 'area_sqft', 'broad_zone', 'place', 'far']]
    austin_parcels = austin_parcels.loc[(austin_parcels['place'] == 'Austin')]
    austin_parcels['percent_undeveloped'] = 1 - austin_parcels['far']
    austin_parcels['developed_sqft'] = np.minimum(austin_parcels['area_sqft'], 
                                                  austin_parcels['far'].multiply(austin_parcels['area_sqft']))

    dallas_path = parcel.get_cached_municipal_parcel_path('dallas')
    dallas_parcels = pd.read_csv(dallas_path)[['dist_to_center', 'area_sqft', 'broad_zone', 'place', 'TOT_MAIN_SF']]
    dallas_parcels = dallas_parcels.loc[(dallas_parcels['place'] == 'Dallas')]
    dallas_parcels['percent_undeveloped'] = 1 - dallas_parcels['TOT_MAIN_SF'].divide(dallas_parcels['area_sqft'])
    dallas_parcels['developed_sqft'] = np.minimum(dallas_parcels['TOT_MAIN_SF'], dallas_parcels['area_sqft'])

    houston_path = parcel.get_cached_municipal_parcel_path('houston')
    houston_parcels = pd.read_csv(houston_path)[['dist_to_center', 'area_sqft', 'broad_zone', 'place', 'BASE_AREA']]
    houston_parcels = houston_parcels.loc[(houston_parcels['place'] == 'Houston')]
    houston_parcels['percent_undeveloped'] = 1 - houston_parcels['BASE_AREA'].divide(houston_parcels['area_sqft'])
    houston_parcels['developed_sqft'] = np.minimum(houston_parcels['BASE_AREA'], houston_parcels['area_sqft'])

    # Group
    all_data = pd.concat([austin_parcels, dallas_parcels, houston_parcels], axis = 0, ignore_index = True)
    all_data['dist_to_center'] = all_data['dist_to_center'].astype(float).apply(np.ceil)
    grouped_data = all_data.groupby(['place', 'dist_to_center', 'broad_zone'])

    # Calculate medians, simple & weighted mean
    medians = grouped_data['percent_undeveloped'].median().reset_index()
    medians['percent_undeveloped'] = 100*medians['percent_undeveloped']
    medians = medians.rename(columns={'percent_undeveloped': 'percent_undeveloped_median'})

    means = grouped_data['percent_undeveloped'].mean().reset_index()
    means['percent_undeveloped'] = 100*means['percent_undeveloped']
    means = means.rename(columns = {'percent_undeveloped': 'percent_undeveloped_simple_mean'})

    total_areas = grouped_data['area_sqft'].sum()
    total_developed_area = grouped_data['developed_sqft'].sum(skipna = True)
    weighted_means = total_developed_area.divide(total_areas).reset_index()
    weighted_means[0] = 100*weighted_means[0]
    weighted_means = weighted_means.rename(columns = {0: 'percent_undeveloped_weighted_average'})

    # Combine and save as csv
    all_results = reduce(lambda left, right: pd.merge(left, right, 
                                                     on=['broad_zone', 'dist_to_center', 'place'], 
                                                     how='outer', 
                                                     sort=False),
                        [weighted_means, means, medians])
    all_results.to_csv('shared_data/calculations/percent_undeveloped.csv')

        # Now graph
    all_results = all_results.loc[all_results['broad_zone'] == 'Single Family']
    all_results = all_results.loc[all_results['dist_to_center'] < 11]
    p = (ggplot(all_results,
                aes(x='dist_to_center', y='percent_undeveloped_simple_mean', fill='place', group='place'))
         + geom_col(position='dodge')
         + labs(title='Undeveloped Area of Single Family Lots in Texas Triangle',
                x='Distance from City Center',
                y='Simple Mean of Undeveloped Percentage of SF Lots'))
    p.save(save_path, width=10, height=8)

def calc_parking_costs(save_path = 'Figures/property_value_histogram.svg'):
    """Calculates average land costs within 1 mile of the city center in Austin, Dallas, Houston. Relies on cached
    municipal parcel data. These are a bit conservative figures because they use lot size instead of base area of the
    actual building."""

    # Read data, get value
    austin_path = parcel.get_cached_municipal_parcel_path('austin')
    austin_parcels = pd.read_csv(austin_path)[['dist_to_center', 'area_sqft', 'place', 'appraised']]
    austin_parcels['value'] = austin_parcels['appraised'].astype(float)
    austin_parcels = austin_parcels.loc[
        (austin_parcels['place'] == 'Austin') & (austin_parcels['dist_to_center'] <= 1)]

    dallas_path = parcel.get_cached_municipal_parcel_path('dallas')
    dallas_parcels = pd.read_csv(dallas_path)[['dist_to_center', 'area_sqft', 'place', 'TOT_VAL']]
    dallas_parcels['value'] = dallas_parcels['TOT_VAL'].astype(float)
    dallas_parcels = dallas_parcels.loc[
        (dallas_parcels['place'] == 'Dallas') & (dallas_parcels['dist_to_center'] <= 1)]

    houston_path = parcel.get_cached_municipal_parcel_path('houston')
    houston_parcels = pd.read_csv(houston_path)[['dist_to_center', 'area_sqft', 'place', 'TOTAL_APPRAISED_VALUE']]
    houston_parcels['value'] = houston_parcels['TOTAL_APPRAISED_VALUE'].astype(float)
    houston_parcels = houston_parcels.loc[
        (houston_parcels['place'] == 'Houston') & (houston_parcels['dist_to_center'] <= 1)]

    # Join, calculate value_per_sqft, group
    all_data = pd.concat([austin_parcels, dallas_parcels, houston_parcels], axis=0)
    all_data['value_per_sqft'] = all_data['value'].divide(all_data['area_sqft'])
    grouped_data = all_data.groupby(['place'])

    # Calculate
    median = grouped_data['value_per_sqft'].median().reset_index()
    median = median.rename(columns={'value_per_sqft': 'value_per_sqft_median'})

    simple_mean = grouped_data['value_per_sqft'].mean().reset_index()
    simple_mean = simple_mean.rename(columns={'value_per_sqft': 'value_per_sqft_simple_mean'})

    weighted_mean = grouped_data['value'].sum().divide(grouped_data['area_sqft'].sum()).reset_index()
    weighted_mean = weighted_mean.rename(columns={0: 'value_per_sqft_mean'})

    all_results = reduce(lambda left, right: pd.merge(left, right, on=['place'], how='outer', sort=False),
                         [weighted_mean, simple_mean, median])
    all_results.to_csv('shared_data/calculations/values_per_sqft.csv')

    # Graph histogram --
    # Get rid of the outliers
    maximum = all_data['value_per_sqft'].quantile(.95)
    all_data.loc[all_data['value_per_sqft'] > maximum, 'value_per_sqft'] = maximum
    all_data = all_data.loc[(all_data['place'].notnull()) & (all_data['value_per_sqft'].notnull())]

    p = (ggplot(all_data, aes(x = 'value_per_sqft', fill = 'place', group = 'place'))
          + geom_histogram()
          + labs(title = 'Property Values within 1 Mile of City Center in Texas Triangle',
                 x = 'Value per Square Foot',
                 y = 'Number of Parcels')
          + theme_bw()
          + facet_wrap('~place'))
    p.save(save_path, width = 10, height = 8)