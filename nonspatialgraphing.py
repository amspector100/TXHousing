import numpy as np
import datetime as dt
import time
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import math
from tqdm import tqdm
import matplotlib
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import helpers
import inputs
import spatial_functions as sf

austin_processed_parcel_data_path = 'data/processed_parcel_data.shp'

# Hist of property values of downzoning and construction vs the city at large for Austin specifically.
def process_parcel_data(zip_intersections = True, save = True):
    """
    Processes austin parcel data. Normalizes property values by land use, computes zip code intersections if zip_intersections = True,  and formats base
    zone properly, and then returns. Will save as a shapefile if save = True.
    :return:
    """

    time0 = time.time()
    print('Starting to read parcel data')
    parcel_data = gpd.read_file(inputs.austin_parcel_path)
    print('Finished reading parcel data, took {} seconds. Now processing base zones.'.format(time.time() - time0))

    # Get basezone in proper format
    def format_basezone(s):
        s = str(s)
        try:
            return s.split('(')[1][:-1]
        except:
            # Usually this means there's no basezone information
            return s
    parcel_data['basezone'] = parcel_data['basezone'].apply(format_basezone)

    # Compute value per area and average value per area for each land use, so we can normalize later
    print('Finished processing base zones, took {} seconds. Now normalizing data conditionally on land use.'.format(time.time() - time0))
    parcel_data.loc[:, 'value_per_area'] = parcel_data['appraised'].divide(parcel_data['shape_area'])
    land_use_dic = {}
    for land_use in parcel_data['land_use'].unique():
        subset = parcel_data.loc[parcel_data['land_use'] == land_use]
        land_use_dic[land_use] = [subset['value_per_area'].mean(), subset['value_per_area'].std()]

    # Normalize - there should be no key areas because we just called parcel_data['land_use'].unique().
    norm_pva = pd.Series(index = parcel_data.index)
    for index, row in parcel_data[['value_per_area', 'land_use']].iterrows():
        norm_pva[index] = (row['value_per_area'] - land_use_dic[row['land_use']][0])/land_use_dic[row['land_use']][1]
    parcel_data['normalized_value_per_area'] = norm_pva


    # Zip intersections
    if zip_intersections:
        print('Finished normalizing, took {} seconds. Computing zip code intersections.'.format(
            time.time() - time0))
        parcel_data = sf.zip_intersect(parcel_data, zips = inputs.austin_zips)


    if save:
        print('Finished processing parcel data, took {} seconds. Now writing to shapefile.'.format(time.time() - time0))
        parcel_data.to_file(austin_processed_parcel_data_path)
        print('Finished writing to shapefile, took {} seconds.'.format(time.time() - time0))
    else:
        print('Finished processing parcel data, took {} seconds.'.format(time.time() - time0))

    return parcel_data


def austin_development_histogram(parcel_data, plot = True):

    time0 = time.time()

    # Get new construction and downzoning subsets
    new_construction = parcel_data.loc[parcel_data['yr_built'] > 2009.0]

    def is_down_developped(old_lu, new_lu):
        """
        Checks whether old land use was a higher residential density than new land use
        """
        dictionary = {160:1, 100:2, 150:3, 210:4, 330:5, 220:6}
        try:
            if dictionary[old_lu] > dictionary[new_lu]:
                return True
            else:
                return False
        except:
            # We will consider it downzoned if it was residential and is no longer
            if old_lu in dictionary and new_lu not in dictionary:
                return True
            else:
                return False

    downzone = parcel_data.loc[parcel_data['previouslu'] != parcel_data['land_use']]
    downzone.loc[:, 'flag'] = False
    for index in downzone.index:
        downzone.at[index, 'flag'] = is_down_developped(downzone.loc[index, 'previouslu'], downzone.loc[index, 'land_use'])
    downzone = downzone.loc[downzone['flag'] == True]

    print('Finished subsetting, took {} seconds.'.format(time.time() - time0))

    if plot:
        print('Now plotting construction and downzone data histogram.')
        # Plot histogram of property values
        bins = np.linspace(-0.9, 3, 40)
        plt.hist(downzone['normalized_value_per_area'].dropna().values, label = 'Downzoned areas', bins = bins, alpha = 0.5, color = 'orange', normed = 1)
        plt.hist(new_construction['normalized_value_per_area'].dropna().values, label = 'New construction', bins = bins, alpha = 0.5, color = 'green', normed = 1)
        plt.hist(parcel_data['normalized_value_per_area'].dropna().values, label = 'All parcel data', bins = bins, alpha = 0.5, color = 'blue', normed = 1)
        plt.legend()
        plt.title('Normalized Property Value Distributions (Conditional on Land Use)')


        # Make y axis prettier by making it percentages instead of standard deviations (I think that's what they were originally)
        def to_percent(y, position):
            # Ignore the (by default) passed in location - scale the default tick locations to percentages
            s = str(100 * y)

            # If latex is on we have to be careful
            if matplotlib.rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                return s + '%'

        formatter = FuncFormatter(to_percent)
        plt.gca().yaxis.set_major_formatter(formatter)

        # Show
        plt.savefig('Figures/Bucket 2/Austin_Downzoning_Hist.png', bbox_inches = 'tight')
        plt.show()

    return downzone, new_construction

if __name__ == '__main__':


    # Testing
    parcel_data = process_parcel_data(zip_intersections=True)
