import numpy as np
import pandas as pd
import geopandas as gpd
import spatial_functions as sf

# Read data from 2000 and 2016 and merge
demo_2000 = pd.read_csv('data/2000-demographics.csv',  low_memory = False)
ec_2000 = pd.read_csv('data/2000-economics.csv', low_memory = False)
dt_2000 = demo_2000.merge(ec_2000, how = 'inner', on = 'Id2')

print(dt_2000.columns[0:5])

pop_2016 = pd.read_csv('data/2016-demographics.csv',  low_memory = False)
housing_2016 = pd.read_csv('data/2016-housing.csv',  low_memory = False)
dt_2016 = housing_2016.merge(pop_2016, how = 'inner', on = 'Id2')

# Reallocate 2016 population to 2000 census tracts. 
# 2000 tracts that were merged in 2010 get proportions of the 2016
# population commensurate with their proportions of the 2010 population.

crosswalk = pd.read_csv('data/tx-census-tract-crosswalk.csv')