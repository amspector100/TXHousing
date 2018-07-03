from inputs import *
from helpers import *
import spatial_functions as sf

# Special setback and minimum lot zone calculations ----------------------------- (Houston) ----------------------------------------------------- Special setback and minimum lot zone calculations
special_setbacks = gpd.read_file(houston_spec_min_setbacks)
special_lots = gpd.read_file(houston_spec_min_lots)

# Find areas
special_setbacks = sf.get_area_in_units(special_setbacks) # Remember, this defaults to miles
special_lots = sf.get_area_in_units(special_lots)
print(special_setbacks['area'].sum())
print(special_lots['area'].sum())

# Find population in the particular areas - B01001e1 is the column which counts population in he block data
block_data = sf.get_block_geodata(['X01_AGE_AND_SEX'], cities = 'Houston')
special_lots = sf.get_all_averages_by_area(block_data, special_lots, fillna = 0) # Defaults to getting pop data
special_setbacks = sf.get_all_averages_by_area(block_data, special_setbacks, fillna = 0)

print(special_lots['B01001e1'].sum())
print(special_setbacks['B01001e1'].sum())


