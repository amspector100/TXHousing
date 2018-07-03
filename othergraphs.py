from inputs import *
from helpers import *
from zipcode_scrape import *
from shapely.ops import cascaded_union
import spatial_functions as sf
import choropleth
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

from plotnine import * # This won't pollute the namespace I don't think because it's exactly like ggplot2 commands

# These booleans determine which graphs to graph in this run of the script. If you want to regenerate every graph,
# just set everything to 'True'.
plot = True
save = True
min_lot_size = False # Min lot size residential areas city
far = False
coverage = False
plot_austin_permits = False
plot_dallas_permits = False
plot_austin_commercial_permits = False
plot_dallas_commercial_permits = False
plot_hds_locations = False
plot_hds_permits = False # This needs some updates/revisions
broad_zone = False
permit_scatter = True


# Read data --------------------------------------------------
dallas_zones = process_zoning_shapefile(dallas_inputs, broaden = True)
austin_zones = process_zoning_shapefile(austin_inputs, broaden = True)
time0 = time.time()

# Minimum lot size graph -------------------------------------
if min_lot_size:

    # Currently excluding agricultural zones ( A(A) in Dallas, RR in Austin) as well as non-residential zones
    dallas_zones['inv_density'] = get_regulation_data(dallas_zones[dallas_inputs.feature], dallas_regulations_path, 'inv_density', fill = 0)
    austin_zones['min_lot'] = get_regulation_data(austin_zones[austin_inputs.feature], austin_regulations_path, 'min_lot', fill = 0)

    print('Min lot size, intersection 1')
    dallas_minlot_radii = sf.polygons_intersect_rings(dallas_zones.loc[(dallas_zones[dallas_inputs.feature] != 'A(A)') &
                                                                       (dallas_zones['broad_zone'] != 'Other')],
                                                      dallas_inputs, factor='inv_density', categorical=False,
                                                                            newproj='epsg:2276', step=1, maximum = 17)

    print('Min lot size, intersection 2')
    austin_minlot_radii = sf.polygons_intersect_rings(austin_zones.loc[(austin_zones[austin_inputs.feature] != 'RR') &
                                                                       (austin_zones['broad_zone'] != 'Other')],
                                                      austin_inputs, factor='min_lot', categorical=False,
                                                                            newproj='epsg:2277', step=1, maximum = 17)

    minlot = pd.concat([dallas_minlot_radii, austin_minlot_radii], axis = 1)
    minlot.columns = ['Dallas', 'Austin']
    minlot.plot(kind='bar', legend = True)
    plt.title('Minimum Lot Size in Austin and Dallas, Excluding Agricultural and Nonresidential Land')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('Average Lot Size (Square Feet)')
    if save:
        plt.savefig('Figures/Bucket 2/Minimum_Lot_Size_Residential_No_Agriculture.png', bbox_inches = 'tight')
    if plot:
        plt.show()

# far graph -------------------------------------

if far:
    dallas_zones['far'] = get_regulation_data(dallas_zones[dallas_inputs.feature], dallas_regulations_path, 'far', fill = None)
    austin_zones['far'] = get_regulation_data(austin_zones[austin_inputs.feature], austin_regulations_path, 'far', fill = None)

    print('Far, intersection 1')
    dallas_far_radii = sf.polygons_intersect_rings(dallas_zones, dallas_inputs, factor='far', categorical=False,
                                                                            newproj='epsg:2276', step=0.5, maximum = 5)

    print('Far, intersection 2')
    austin_far_radii = sf.polygons_intersect_rings(austin_zones, austin_inputs, factor='far', categorical=False,
                                                                            newproj='epsg:2277', step=0.5, maximum = 5)

    fardata = pd.concat([dallas_far_radii, austin_far_radii], axis = 1)
    fardata.columns = ['Dallas', 'Austin']
    fardata.plot(kind='bar', legend = True)
    plt.title('Floor to Area Ratio in Austin and Dallas')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('Floor to Area Ratio')
    if save:
        plt.savefig('Figures/Bucket 2/FAR.png', bbox_inches = 'tight')
    if plot:
        plt.show()

# coverage --------------------------------------------------------
if coverage:
    dallas_zones['coverage'] = get_regulation_data(dallas_zones[dallas_inputs.feature], dallas_regulations_path, 'coverage', fill = 100)
    austin_zones['max_build_cov'] = get_regulation_data(austin_zones[austin_inputs.feature], austin_regulations_path, 'max_build_cov', fill = 100)

    print(dallas_zones['coverage'])

    print('cov, intersection 1')
    cov_dallas = sf.polygons_intersect_rings(dallas_zones, dallas_inputs, factor='coverage', categorical=False,
                                                                            newproj='epsg:2276', step=1, maximum = 19)

    print('cov, intersection 2')
    cov_austin = sf.polygons_intersect_rings(austin_zones, austin_inputs, factor='max_build_cov', categorical=False,
                                                                            newproj='epsg:2277', step=1, maximum = 19)

    coverage_data = pd.concat([cov_dallas, cov_austin], axis = 1)
    coverage_data.columns = ['Dallas', 'Austin']
    coverage_data.plot(kind='bar', legend = True)
    plt.title('Maximum Building Coverage in Austin and Dallas')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('Max Building Coverage (%)')
    if save:
        plt.savefig('Figures/Bucket 2/BldgCoverage.png', bbox_inches = 'tight')
    if plot:
        plt.show()

# Permits for single/multifamily building
if plot_dallas_permits:

    dallas_permit_data = get_corrected_dallas_permit_data() #Get dallas permit data from the saved file with CORRECTED geocodes and the original columns.

    dallas_permit_rings = sf.points_intersect_rings(dallas_permit_data, dallas_inputs, factor = 'Permit Type', step = 1, categorical = True, by = 'median', per_square_mile = False)
    dallas_permit_rings.plot(kind = 'bar', legend = True)
    plt.title('Residential Permits Issued, Dallas, TX')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('New Construction Permits Issued, 2011-2016')
    if save:
        plt.savefig('Figures/Bucket 2/DallasResPermits.png', bbox_inches = 'tight')
    if plot:
        plt.show()

if plot_austin_permits:

    austin_permit_data = process_austin_permit_data(searchfor=['101 single family houses',
                                          '103 two family bldgs',
                                          '104 three & four family bldgs',
                                          '105 five or more family bldgs'],
                               earliest=2013)


    def map_permitclass(text):
        if '101' in text:
            return 'Single Family Permit'
        else:
            return 'Mutifamily Permit'

    austin_permit_data['PermitClass'] = austin_permit_data['PermitClass'].apply(map_permitclass)

    austin_permit_rings = sf.points_intersect_rings(austin_permit_data, austin_inputs, factor = 'PermitClass', step = 1, categorical = True, by = 'median')
    austin_permit_rings.plot(kind = 'bar', legend = True)
    plt.title('Residential Permits Issued per Square Mile, Austin, TX')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('New Construction Permits Issued Per Square Mile, 2013-2018')
    if save:
        plt.savefig('Figures/Bucket 2/AustinResPermitsPSM.png', bbox_inches = 'tight')
    if plot:
        plt.show()

if plot_austin_commercial_permits:

    pass

if plot_hds_locations or plot_hds_permits:

    # Make zones valid
    austin_zones = sf.process_geometry(austin_zones)

    # Historic Districts (hds)

    # Get national hds and transform
    tx_hd_path = "data/Zoning Shapefiles/NationalRegisterPY_shp/NationalRegisterPY.shp"
    tx_hd_data = gpd.read_file(tx_hd_path)
    tx_hd_data = tx_hd_data.to_crs({'init':'epsg:4326'})

    # Get austin hds and others
    austin_signature = '-HD'
    austin_local_hds = austin_zones.loc[[austin_signature in text for text in austin_zones[austin_inputs.feature]]]
    austin_nat_hds = tx_hd_data.loc[tx_hd_data['CITY'] == 'Austin']

    # Get Dallas hds and others
    dallas_zones = dallas_zones.to_crs({'init':'epsg:4326'})
    dallas_cds = dallas_zones.loc[dallas_zones['LONG_ZONE_'].apply(lambda x: x[0:2]) == 'CD']
    dallas_nat_hds = tx_hd_data.loc[tx_hd_data['CITY'] == 'Dallas']

    # Now, just plot distance from city center ---------- - - - - - -
    if plot_hds_locations:
        step = 0.5
        maximum = 5
        austin_local_data =  sf.polygons_intersect_rings(austin_local_hds, austin_inputs, factor=None, newproj='epsg:2276', step=step, maximum=maximum)
        austin_local_data.name = 'Austin Local Historic Districts'
        austin_nat_data = sf.polygons_intersect_rings(austin_nat_hds, austin_inputs, factor=None, newproj='epsg:2276', step=step, maximum=maximum)
        austin_nat_data.name = 'Austin National Historic Districts'

        dallas_cds_data = sf.polygons_intersect_rings(dallas_cds, dallas_inputs, factor=None, newproj='epsg:2276', step=step, maximum=maximum)
        dallas_cds_data.name = 'Dallas Conservation Districts'
        dallas_nat_data = sf.polygons_intersect_rings(dallas_nat_hds, dallas_inputs, factor=None, newproj='epsg:2276', step=step, maximum=maximum)
        dallas_nat_data.name = 'Dallas Historic Districts'


        all_loc_data = pd.concat([austin_local_data, austin_nat_data, dallas_cds_data, dallas_nat_data], axis = 1)
        all_loc_data.plot(kind = 'bar', color = ['purple', 'blue', 'red', 'orange'], legend = True)
        plt.title('Locations of Austin and Dallas Historic Zones in relation to the City Center')
        plt.xlabel('Distance from the City Center, Miles')
        plt.ylabel('Percent of Land Covered by District Type')
        if save:
            plt.savefig('Figures/Bucket 2/HDLocations.png', bbox_inches = 'tight')
        if plot:
            plt.show()

    # Now plot number of demolition permits in the areas ------------------ - -- - - - -
    if plot_hds_permits:
        # Get Renovation/Demolition (not construction) points/data.
        dallas_points = process_dallas_permit_data(permit_types = ['Demolition Permit Commercial', 'Demolition Permit SFD/Duplex']) #dallas_renovation_types
        austin_points = process_austin_permit_data(searchfor = ['demolition'], earliest = 2013)
        austin_points = austin_points.loc[[not bool for bool in austin_points['Description'].str.contains('interior')]] # Don't care about interior remodels
        print('Length of austin_points is {}'.format(len(austin_points)))

        radius = 2

        # Get residential cores. Need to simplify slightly to save a massive amount of time
        dallas_core = sf.get_urban_core(dallas_inputs, radius)
        print('Taking complex Dallas urban core intersection, time is {}'.format(time.time()-time0))
        dallas_res_core = dallas_core.intersection(dallas_zones.loc[dallas_zones['broad_zone'] != 'Other', 'geometry'].unary_union)
        dallas_res_core = gpd.GeoDataFrame(geometry = dallas_res_core)
        dallas_res_core.loc[:, 'geometry']  = dallas_res_core['geometry'].simplify(tolerance = 0.001)

        austin_core = sf.get_urban_core(austin_inputs, radius)
        print('Taking complex Austin urban core intersection, time is {}'.format(time.time() - time0))
        austin_res_core = austin_core.intersection(austin_zones.loc[austin_zones['broad_zone'] != 'Other', 'geometry'].unary_union)
        austin_res_core = gpd.GeoDataFrame(geometry = austin_res_core)
        austin_res_core.loc[:, 'geometry']  = austin_res_core['geometry'].simplify(tolerance = 0.001)


        # Get everything else
        def points_over_area(points, polygons, points_geometry_column = 'geometry', poly_geometry_column = 'geometry', normalize_by_area = True):
            """
            :param points: GeoDataframe of points, should be in lat long
            :param polygons: Geodataframe of polygons
            :param newproj: Transform data into this new projection to calculate areas more effectively. Units must be in feet.
            :param: normalize_by_area: If true, will divide the number of points by the area of the location (in square miles).
            :return: GeoDataFrame of polygons
            """

            # Process points and polys to make them valid
            polygons = sf.process_geometry(polygons, drop_multipolygons = False)
            polygons.reset_index(drop=True, inplace=True)
            polygons.index = [str(ind) for ind in polygons.index]

            points = points.loc[points[points_geometry_column].is_valid] #Ignore invalid points (i.e. with 'na's)
            points.reset_index(drop=True, inplace=True)
            points.index = [str(ind) for ind in points.index]

            counter = 0

            # Intersections (efficiently)
            spatial_index = points.sindex
            for poly in polygons[poly_geometry_column]:
                possible_intersections_index = list(spatial_index.intersection(poly.bounds))
                possible_intersections = points.iloc[possible_intersections_index]
                precise_matches = possible_intersections[points_geometry_column].intersects(poly)
                counter += sum(precise_matches.tolist())

            # Return # of points divided by area
            area = sf.get_area_in_units(polygons)['area'].sum()

            if normalize_by_area:
                return counter/area # Would be nice to get this area in miles. Oh well.
            else:
                return counter

        # Label names
        local_name = 'Conservation District (Dallas) or Local Historic Distrct (Austin)'
        nat_name = 'National Historic District'
        all_name = 'Urban Core ({} Miles from City Center)'.format(radius)

        # Initialize and fill results
        result = pd.DataFrame(index = ['Dallas', 'Austin'], columns = [local_name, nat_name, all_name])


        # Get block data
        block_data = sf.get_block_geodata(['X01_AGE_AND_SEX'], cities = ['Austin', 'Dallas'])


        print('Filling hd result for Dallas, time is {}'.format(time.time() - time0))
        for name, data in zip([local_name, nat_name], [dallas_cds, dallas_nat_hds]):
            result.at['Dallas', name] = points_over_area(dallas_points, data, normalize_by_area = False)/(sf.get_all_averages_by_area(block_data['Dallas'], data, drop_multipolygons = False, fillna = 0, account_for_water = True)['B01001e1'].sum())

        # Do this separately to save time in getting num_points
        num_points = points_over_area(dallas_points,
                                     dallas_core,
                                     normalize_by_area=False)
        pop = sf.get_all_averages_by_area(block_data['Dallas'],
                                           dallas_res_core,
                                           drop_multipolygons = False,
                                           fillna = 0,
                                           account_for_water = True)['B01001e1'].sum()
        result.at['Dallas', all_name] = num_points/pop

        print('Filling result for Austin, time is {}'.format(time.time() - time0))
        for name, data in zip([local_name, nat_name], [austin_local_hds, austin_nat_hds]):
            # Get number of points per area and then also divide by population
            result.at['Austin', name] = points_over_area(austin_points, data, normalize_by_area = False)/(sf.get_all_averages_by_area(block_data['Austin'], data, drop_multipolygons = False, fillna = 0, account_for_water = True)['B01001e1'].sum())
        # Do this separately to save time in getting num_points
        num_points = points_over_area(austin_points,
                                     austin_core,
                                     normalize_by_area=False)
        pop = sf.get_all_averages_by_area(block_data['Austin'],
                                           austin_res_core,
                                           drop_multipolygons = False,
                                           fillna = 0,
                                           account_for_water = True)['B01001e1'].sum()
        result.at['Austin', all_name] = num_points/pop


        result = 100*result.divide(result[all_name], axis = 0)
        result.plot(kind = 'bar', legend = True)
        plt.title('Demolition Permit Frequencies in the Urban Core and in Historic Districts, Austin and Dallas')
        plt.xlabel('District Type')
        plt.ylabel('Demolition Permits per Population, Scaled to Urban Core (%)')
        if save:
            plt.savefig('Figures/Bucket 2/HDDemolitionPermits.png', bbox_inches = 'tight')
        if plot:
            plt.show()

if broad_zone:

    # This currently excludes agricultural and nonresidential land.
    maximum = 17

    dallas_zones = process_zoning_shapefile(dallas_inputs, broaden=True)
    austin_zones = process_zoning_shapefile(austin_inputs, broaden=True)

    dallas_zone_rings = sf.polygons_intersect_rings(dallas_zones.loc[(dallas_zones[dallas_inputs.feature] != 'A(A)') &
                                                                       (dallas_zones['broad_zone'] != 'Other')],
                                                      dallas_inputs, factor='broad_zone', categorical=True,
                                                                            newproj='epsg:2276', step=1, maximum = maximum)

    dallas_zone_rings.columns = [str(col) + ' Dallas' for col in dallas_zone_rings.columns]




    austin_zone_rings = sf.polygons_intersect_rings(austin_zones.loc[(austin_zones[austin_inputs.feature] != 'RR') &
                                                                       (austin_zones['broad_zone'] != "Other")],
                                                      austin_inputs, factor='broad_zone', categorical=True,
                                                                            newproj='epsg:2277', step=1, maximum = maximum)
    austin_zone_rings.columns = [str(col) + ' Austin' for col in austin_zone_rings.columns]



    zone_rings = pd.concat([dallas_zone_rings, austin_zone_rings], axis = 1)
    zone_rings[['Single Family Austin', 'Single Family Dallas']].plot(kind='bar', legend = True, color = ['cornflowerblue', '#8FBC8F'])
    plt.title('Single Family Zoning in Austin and Dallas, Excluding Agricultural and Nonresidential Land')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('Percentage of Residential Land Zoned')
    if save:
        plt.savefig('Figures/Bucket 2/SFZoningRings.png', bbox_inches = 'tight')
    if plot:
        plt.show()

    zone_rings[['Multifamily Austin', 'Multifamily Dallas']].plot(kind='bar', legend = True, color = ['cornflowerblue', '#8FBC8F'])
    plt.title('Multifamily Zoning in Austin and Dallas, Excluding Agricultural and Nonresidential Land')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('Percentage of Residential Land Zoned')
    if save:
        plt.savefig('Figures/Bucket 2/MFZoningRings.png', bbox_inches = 'tight')
    if plot:
        plt.show()

if permit_scatter:

    # Get permit data
    dallas_permits = get_corrected_dallas_permit_data(path = dpm_save_path)
    austin_permits = process_austin_permit_data(searchfor=['101 single family houses',
                                                           '103 two family bldgs',
                                                           '104 three & four family bldgs',
                                                           '105 five or more family bldgs'],
                                                earliest=2013,
                                                permittypedesc='Building Permit',
                                                workclass='New')

    # Now get zip geodata - this is for the scatter of med housing price
    zip_geodata = get_zip_boundaries()
    austin_zip_geodata = zip_geodata.loc[[z for z in austin_zips if z in zip_geodata.index]]
    dallas_zip_geodata = zip_geodata.loc[[z for z in dallas_zips if z in zip_geodata.index]]
    # Loop through and add pricing data
    for cityname, citydata in zip(['Autin, TX', 'Dallas, TX'], [austin_zip_geodata, dallas_zip_geodata]):
        citydata = add_demand_data(zip_geodata=citydata, demand_input = realtor_avg_sf_price, city = cityname, feature_name = 'sf_avg_listing')
        citydata = add_demand_data(zip_geodata=citydata, demand_input = realtor_avg_cth_price, city = cityname, feature_name = 'mf_avg_listing')

    # Now find which points are in which zipcodes
    for sfpermits, mfpermits, data in zip([austin_permits, dallas_permits], [austin_zip_geodata, dallas_zip_geodata]):
        spatial_index = permits.sindex
        # Initialize result column
        permits['zipcode'] = None
        for zipcode, polygon in tqdm(zip(data.index, data['geometry'])):
            possible_matches_index = list(spatial_index.intersection(polygon.bounds))
            possible_matches = permits.iloc[possible_matches_index]
            permits.loc[possible_matches['geometry'].intersects(polygon), 'zipcode'] = zipcode

    print(austin_permit_data)


    austin_sf = austin_permits.loc[austin_permits['PermitClass'].str.contains('101 single')]
    austin_mf = austin_permits.loc[[not bool for bool in austin_permits['PermitClass'].str.contains('101 single')]]

    dallas_sf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Single')]
    dallas_mf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Multi')]

