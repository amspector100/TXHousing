from inputs import *
from helpers import *
import spatial_functions as sf
from tqdm import tqdm

special_setbacks_and_lots = False
calc_percent_residential = False
calc_landmarks_in_houston = True
calc_parking_cost = False

# Special setback and minimum lot zone calculations ----------------------------- (Houston) ----------------------------------------------------- Special setback and minimum lot zone calculations
if special_setbacks_and_lots:

    print('Running some numbers on special setbacks and lots')
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


# Calculate percent of land that is residential ------------------------------ (Austin/Dallas) ---------------------------------------------------- Calculate percent of land that is residential

if calc_percent_residential:
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
    results.to_csv('data/bg_percent_residential.csv')

# Calculate population of nat hist districts in Houston as well as
if calc_landmarks_in_houston:

    # Get points and polygons
    houston_landmarks = gpd.read_file(houston_historic_landmarks_path).to_crs({'init':'epsg:4326'})
    tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
    tx_hd_data = gpd.read_file(tx_hd_path)
    houston_nat_districts = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston'].to_crs({'init':'epsg:4326'})
    houston_block_data = sf.get_block_geodata(['X01_AGE_AND_SEX'], cities = 'Houston')
    fig, ax = plt.subplots()
    houston_block_data.plot(ax = ax, column = 'B01001e1', edgecolor = 'black', legend = True)
    houston_nat_districts.plot(ax = ax, color = 'blue', alpha = 0.5)
    plt.show()

    # Calculate number of points in polygons
    num_points = sf.points_over_area(houston_landmarks, houston_nat_districts, normalize_by_area = False)
    population = sf.get_all_averages_by_area(houston_block_data, houston_nat_districts, fillna=0,
                                             account_method='water')['B01001e1'].sum()

    print('Number of Historic Landmarks in NRHDS is {}, Population of NRHDs is {}'.format(num_points, population))

# Calculate parking cost
if calc_parking_cost:

    austin_data = gpd.read_file(austin_parcel_path)
    austin_data