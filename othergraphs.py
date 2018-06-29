from inputs import *
from helpers import *
from zipcode_scrape import *
import spatial_functions as sf
import choropleth

plot = True
save = True
min_lot_size = False
far = False
coverage = False
plot_austin_permits = False
plot_dallas_permits = False
plot_austin_commercial_permits = False
plot_dallas_commercial_permits = False
plot_hds = True

broad_zone = False

# Read data --------------------------------------------------
dallas_zones = process_zoning_shapefile(dallas_inputs, broaden = True)
austin_zones = process_zoning_shapefile(austin_inputs, broaden = True)





# Minimum lot size graph -------------------------------------
if min_lot_size:

    # Currently excluding agricultural zones ( A(A) in Dallas, RR in Austin)

    dallas_zones['inv_density'] = get_regulation_data(dallas_zones[dallas_inputs.feature], dallas_regulations_path, 'inv_density', fill = 0)
    austin_zones['min_lot'] = get_regulation_data(austin_zones[austin_inputs.feature], austin_regulations_path, 'min_lot', fill = 0)

    print('Min lot size, intersection 1')
    dallas_minlot_radii = sf.polygons_intersect_rings(dallas_zones.loc[dallas_zones[dallas_inputs.feature] != 'A(A)'],
                                                      dallas_inputs, factor='inv_density', categorical=False,
                                                                            newproj='epsg:2276', step=1, maximum = 19)

    print('Min lot size, intersection 2')
    austin_minlot_radii = sf.polygons_intersect_rings(austin_zones.loc[austin_zones[austin_inputs.feature] != 'RR'],
                                                      austin_inputs, factor='min_lot', categorical=False,
                                                                            newproj='epsg:2277', step=1, maximum = 19)

    minlot = pd.concat([dallas_minlot_radii, austin_minlot_radii], axis = 1)
    minlot.columns = ['Dallas', 'Austin']
    minlot.plot(kind='bar', legend = True)
    plt.title('Minimum Lot Size in Austin and Dallas, Excluding Agricultural Land')
    plt.xlabel('Distance from the city center (Miles)')
    plt.ylabel('Average Lot Size (Square Feet)')
    if save:
        plt.savefig('Figures/Bucket 2/Minimum_Lot_Size_No_Agriculture.png', bbox_inches = 'tight')
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

if plot_hds:

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

    # Get RENOVATION (not construction) points/data.
    dallas_points = process_dallas_permit_data(permit_types = ['Demolition Permit Commercial', 'Demolition Permit SFD/Duplex'])
    austin_points = process_austin_permit_data(searchfor = ['demolition'])

    radius = 5

    # Get only residential areas
    austin_core = sf.get_urban_core(austin_inputs, radius)
    dallas_core = sf.get_urban_core(dallas_inputs, radius)

    # Get everything else
    def points_over_area(points, polygons, points_geometry_column = 'geometry', poly_geometry_column = 'geometry', newproj = 'epsg:2277', normalize_by_area = True):
        """
        :param points: GeoDataframe of points, should be in lat long
        :param polygons: Geodataframe of polygons
        :param newproj: Transform data into this new projection to calculate areas more effectively. Units must be in feet.
        :param: normalize_by_area: If true, will divide the number of points by the area of the location (in square miles).
        :return: GeoDataFrame of polygons
        """

        # Process points and polys to make them valid
        polygons = sf.process_geometry(polygons)
        polygons.reset_index(drop=True, inplace=True)
        polygons.index = [str(ind) for ind in polygons.index]

        points = points.loc[points[points_geometry_column].is_valid] #Ignore invalid points (i.e. with 'na's)
        polygons.reset_index(drop=True, inplace=True)
        polygons.index = [str(ind) for ind in polygons.index]

        counter = 0

        # Intersections (efficiently)
        spatial_index = points.sindex
        for poly in polygons[poly_geometry_column]:
            possible_intersections_index = list(spatial_index.intersection(poly.bounds))
            possible_intersections = points.iloc[possible_intersections_index]
            precise_matches = possible_intersections[points_geometry_column].intersects(poly)
            counter += sum(precise_matches.tolist())

        # Return # of points divided by area
        polygons = polygons.to_crs({'init':newproj})
        total_area = sum(polygons[poly_geometry_column].area)*(3.58701*10**(-8)) # Square feet to square miles
        if normalize_by_area:
            return counter/total_area # Would be nice to get this area in miles. Oh well.
        else:
            return counter

    # Label names
    local_name = 'Conservation District or Local Historic Distrct'
    nat_name = 'National Historic District'
    all_name = 'Urban Core ({} Miles from City Center)'.format(radius)

    # Initialize and fill results
    result = pd.DataFrame(index = ['Dallas', 'Austin'], columns = [local_name, nat_name, all_name])


    # Get block data
    block_data = sf.get_block_geodata(['X01_AGE_AND_SEX'], cities = ['Austin', 'Dallas'])


    print('Filling result for Austin')
    for name, data in zip(result.columns, [austin_local_hds, austin_nat_hds, austin_core]):
        # Get number of points per area and then also divide by population
        result.at['Austin', name] = points_over_area(austin_points, data, normalize_by_area = False)/(sf.get_all_averages_by_area(block_data['Austin'], data, fillna = 0)['B01001e1'].sum())

    print('Filling hd result for Dallas')
    for name, data in zip(result.columns, [dallas_cds, dallas_nat_hds, dallas_core]):
        result.at['Dallas', name] = points_over_area(dallas_points, data, normalize_by_area = False)/(sf.get_all_averages_by_area(block_data['Dallas'], data, fillna = 0)['B01001e1'].sum())

    result = result.divide(result[all_name], axis = 0)

    result.plot(kind = 'bar', legend = True)
    plt.show()





