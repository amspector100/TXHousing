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
from inputs import austin_processed_parcel_data_path
import spatial_functions as sf

# THis file will eventually be depreciated - it just includes (now unimportant austin and dallas develpment hists.

def austin_development_histogram(parcel_data, plot = True, save = True):

    time0 = time.time()

    # Get new construction and downzoning subsets
    new_construction = parcel_data.loc[parcel_data['yr_built'] > 2011.0] # This is the year

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

    flags = pd.Series(index = downzone.index)
    for index in flags.index:
        flags[index] = is_down_developped(downzone.loc[index, 'previouslu'], downzone.loc[index, 'land_use'])
    downzone = downzone.loc[flags.values]
    print('Finished subsetting, took {} seconds.'.format(time.time() - time0))

    if plot:
        print('Now plotting construction and downzone data histogram.')
        # Plot histogram of property values
        bins = np.linspace(-0.9, 3, 40)
        plt.hist(downzone['normalized_value_per_area'].dropna().values, label = 'Downzoned areas', bins = bins, alpha = 0.5, color = 'orange', normed = 1)
        plt.hist(new_construction['normalized_value_per_area'].dropna().values, label = 'New construction', bins = bins, alpha = 0.5, color = 'green', normed = 1)
        plt.hist(parcel_data['normalized_value_per_area'].dropna().values, label = 'All parcel data', bins = bins, alpha = 0.5, color = 'blue', normed = 1)
        plt.legend()
        plt.title('Austin Normalized Property Value Distributions (Conditional on Land Use and Shape Area)')


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

    if save:
        print('Now saving data.')
        new_construction.to_file(inputs.austin_construction_path)
        downzone.to_file(inputs.austin_downzone_path)

    return downzone, new_construction

def dallas_development_histogram(parcel_data, plot = True, save = True):

    time0 = time.time()

    # Get new construction and downzoning subsets
    construction = parcel_data.loc[(parcel_data['year_blt_2013'] > 2008.0) & (parcel_data['res_com_2013'] == 'R')] # Last 5 years of res construction

    def is_down_developped(old_sp, new_sp):
        """
        Checks whether old sptbcode was a higher residential density than new sptbcode
        """
        dictionary = {'M31':0, 'M32':0, 'A20':0, 'A11':1, 'A12':2, 'A13':3, 'B12':4, 'B11':5,
                      'C11':-1, 'C12':-1, 'C13':-1, 'C14':-1} # -1 represents vacant lots which haven't been downzoned
        try:
            if dictionary[old_sp] > dictionary[new_sp] and dictionary[new_sp] != -1:
                return True
            else:
                return False
        except:
            # We will consider it downzoned if it was residential and is no longer
            if old_sp in dictionary and new_sp not in dictionary:
                return True
            else:
                return False

    downdev = parcel_data.loc[parcel_data['sptbcode_2013'] != parcel_data['sptbcode_2016']]
    flags = pd.Series(index = downdev.index)
    for index in flags.index:
        flags[index] = is_down_developped(downdev.loc[index, 'sptbcode_2013'], downdev.loc[index, 'sptbcode_2016'])
    downdev = downdev.loc[flags.values]

    print('Finished subsetting, took {} seconds.'.format(time.time() - time0))

    if plot:
        print('Now plotting construction and downzone data histogram.')
        # Plot histogram of property values
        bins = np.linspace(-0.9, 3, 40)
        plt.hist(downdev['normalized_value_per_area'].dropna().values, label = 'Downzoned areas', bins = bins, alpha = 0.5, color = 'orange', normed = 1)
        plt.hist(construction['normalized_value_per_area'].dropna().values, label = 'New construction', bins = bins, alpha = 0.5, color = 'green', normed = 1)
        plt.hist(parcel_data['normalized_value_per_area'].dropna().values, label = 'All parcel data', bins = bins, alpha = 0.5, color = 'blue', normed = 1)
        plt.legend()
        plt.title('Dallas Normalized Property Value Distributions (Conditional on SPTBCODE and Shape Area)')


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
        plt.savefig('Figures/Bucket 2/Dallas_Downzoning_Hist.png', bbox_inches = 'tight')
        plt.show()

    if save:
        print('Now saving data.')
        downdev.to_file(inputs.dallas_downdev_path)
        construction.to_file(inputs.dallas_construction_path)

    return downdev, construction





if __name__ == '__main__':


    # Testing
    #parcel_data = process_austin_parcel_data(zip_intersections=True)
    #parcel_data = helpers.process_dallas_parcel_data()
    #dallas_development_histogram(parcel_data, plot = True, save = True)
    parcel_data = process_austin_parcel_data(zip_intersections = False, save = False)
    austin_development_histogram(parcel_data, plot = True, save = True)
