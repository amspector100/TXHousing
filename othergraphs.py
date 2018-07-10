from inputs import *
from helpers import *
from zipcode_scrape import *
from shapely.ops import cascaded_union
import spatial_functions as sf
import choropleth
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import copy

from plotnine import * # This won't pollute the namespace I don't think because it's exactly like ggplot2 commands

import folium
from folium import FeatureGroup
from BindColorMap import *

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
permit_scatter = False
income_histogram = False
construction_heatmap = True

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
        # Update: it's pretty hard to get exterior renovation data for Houston so we're going to focus mostly on
        # demolition permits. This is what we the paper is more interested in anyway.

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

            points = sf.process_points(points, geometry_column = points_geometry_column)

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

    # Get permit data - Dallas
    print('Processing Dallas permit data')
    dallas_permits = get_corrected_dallas_permit_data(path = dpm_save_path)
    dallas_permits = sf.process_points(dallas_permits)
    dallas_sf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Single')]
    dallas_mf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Multi')]

    # Get permit data - Austin
    austin_permits = process_austin_permit_data(searchfor=['101 single family houses',
                                                           '103 two family bldgs',
                                                           '104 three & four family bldgs',
                                                           '105 five or more family bldgs'],
                                                earliest=2013,
                                                permittypedesc='Building Permit',
                                                workclass='New')
    austin_permits = sf.process_points(austin_permits)
    austin_sf = austin_permits.loc[austin_permits['PermitClass'].str.contains('101 single')]
    austin_mf = austin_permits.loc[[not bool for bool in austin_permits['PermitClass'].str.contains('101 single')]]

    # Get permit data - Houston
    # Construction permit data
    print('Processing Houston permit data')
    houston_permit_data = process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE',
                                                                   'NEW AP', 'NEW HI-'],
                                                      searchin = ['PROJ_DESC'],
                                                      earliest = 2013, latest = None)

    # Subset to only include approved permits and nonempty geometries
    houston_permit_data = houston_permit_data.loc[houston_permit_data['Approval'] == 1.0]
    houston_permit_data = sf.process_points(houston_permit_data)
    houston_sf = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW S.F.',
                                                                                                 'NEW SF',
                                                                                                 'NEW TOWNHOUSE',
                                                                                                 'NEW SINGLE']))]
    houston_mf = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW AP',
                                                                                                 'NEW HI-']))]



    # Now get zip geodata - this is for the scatter of med housing price
    zip_geodata = get_zip_boundaries()
    austin_zip_geodata = zip_geodata.loc[[z for z in austin_zips if z in zip_geodata.index]]
    dallas_zip_geodata = zip_geodata.loc[[z for z in dallas_zips if z in zip_geodata.index]]
    houston_zip_geodata = zip_geodata.loc[[z for z in houston_zips if z in zip_geodata.index]]

    # Loop through and add pricing data
    for cityname, citydata in zip(['Austin, TX', 'Dallas, TX', 'Houston, TX'], [austin_zip_geodata, dallas_zip_geodata, houston_zip_geodata]):
        citydata = add_demand_data(zip_geodata=citydata, demand_input = realtor_avg_sf_price, city = cityname, feature_name = 'sfprice')
        citydata = add_demand_data(zip_geodata=citydata, demand_input = realtor_avg_cth_price, city = cityname, feature_name = 'mfprice')
        citydata = add_demand_data(zip_geodata=citydata, demand_input = realtor_hotness_cbsa, city = cityname, feature_name = 'vpp')

    # Now find which points are in which zipcodes
    for sfpermits, mfpermits, data in zip([austin_sf, dallas_sf], [austin_mf, dallas_mf], [austin_zip_geodata, dallas_zip_geodata]):

        # Spatial indexes
        sf_spatial_index = sfpermits.sindex
        mf_spatial_index = mfpermits.sindex

        # Run through
        def get_sf_pts(poly):
            return sf.points_intersect_polygon(sfpermits, poly, sf_spatial_index)
        def get_mf_pts(poly):
            return sf.points_intersect_polygon(mfpermits, poly, mf_spatial_index)

        # I have reality checked this
        data.loc[:, 'SF'] = data['geometry'].apply(get_sf_pts)
        data.loc[:, 'MF'] = data['geometry'].apply(get_mf_pts)

    # Houston already has zip code information listed so we'll just get it
    def num_sf_permits_houston(zipcode):
        return len(houston_sf.loc[houston_sf['Zipcode'] == zipcode])
    def num_mf_permits_houston(zipcode):
        return len(houston_mf.loc[houston_mf['Zipcode'] == zipcode])

    houston_zip_geodata.loc[:, 'SF'] = houston_zip_geodata['ZCTA5CE10'].apply(num_sf_permits_houston)
    houston_zip_geodata.loc[:, 'MF'] = houston_zip_geodata['ZCTA5CE10'].apply(num_mf_permits_houston)


    # Add city information
    houston_zip_geodata['City'] = 'Houston'
    austin_zip_geodata['City'] = 'Austin'
    dallas_zip_geodata['City'] = 'Dallas'

    all_zip_geodata = austin_zip_geodata.append(dallas_zip_geodata).append(houston_zip_geodata)
    y_max = 1000

    #all_zip_geodata = all_zip_geodata.melt(value_vars = ['SF', 'MF'], var_name = 'Housing Type', value_name = 'Number of Permits')
    sfplot = (ggplot(all_zip_geodata, aes(x = 'sfprice', y = 'SF', color = 'City')) #.loc[all_zip_geodata['SF'] <= outlier_threshhold]
              + geom_point()
              + stat_smooth(method = 'lowess', span = 0.5)
              + facet_wrap('City', scales = 'fixed')
              + scale_x_continuous(limits = (0, 2200000), breaks = (0, 1000000, 2000000))
              + scale_y_continuous(limits = (0, y_max))
              + labs(title = 'New Single Family Construction Permits against SF Housing Price',
                     caption = 'Permit Data from the Cities of Austin, Dallas, and Houston /n, Pricing data from Realtor')
              + xlab('Median Single Family Home Price (by Zip Code)')
              + ylab('Number of SF Construction Permits Issued In Past 5 Years'))

    sfplot_variant = (ggplot(all_zip_geodata, aes(x = 'vpp', y = 'SF', color = 'City'))
                      + geom_point()
                      + stat_smooth(method = 'lowess', span = 0.5)
                      #+ scale_x_continuous(limits = (0, 2200000), breaks = (0, 1000000, 2000000))
                      + labs(title = 'New Single Family Construction Permits against Realtor Hotness Index')
                      + xlab('Realtor Hotness Index (by Zip Code)')
                      + ylab('Number of SF Construction Permits Issued in Past 5 Years')
                      + facet_wrap('City', scales='fixed'))

    mfplot = (ggplot(all_zip_geodata, aes(x = 'mfprice', y = 'MF', color = 'City'))
          + geom_point()
          + stat_smooth(method = 'lowess', span = 0.7)
          + facet_wrap('City', scales='free'))

    width = 10
    height = 5
    sfplot.save(filename='Figures/Bucket 2/construction_scatter_sf.svg', width=width, height=height, bbox_inches = 'tight')
    mfplot.save(filename='Figures/Bucket 2/construction_scatter_mf.svg', width=width, height=height, bbox_inches = 'tight')
    sfplot_variant.save(filename='Figures/Bucket 2/construction_scatter_sf_variant.svg', width=width, height=height, bbox_inches='tight')

    # Make a map of the outliers for Connor. Start by plotting the outlier zip codes. ---------------------------------------
    print('Saved files, now working on folium map of permit outliers')

    outlier_threshhold = 0.9 # Quantile
    basemap = folium.Map([austin_inputs.lat, austin_inputs.long], zoom_start = 7)
    for citydata, cityname in zip([austin_zip_geodata, dallas_zip_geodata, houston_zip_geodata], ['Austin', 'Dallas', 'Houston']):
        for factor, mid_color, end_color in zip(['SF', 'MF'], ['blue', 'red'], ['navy', '#660000']):
            min_value = citydata[factor].quantile(outlier_threshhold)
            print(min_value)
            filtered_data = citydata.loc[citydata[factor] >= min_value]
            print(filtered_data)
            layer_name = 'Zip Codes with Top {} % of {} Construction Permits in {}'.format(round(100*(1 - outlier_threshhold)), factor, cityname)
            scale_name = 'Number of New {} Construction Permits in Top {} % of Zipcodes in {}'.format(factor, round(100*(1-outlier_threshhold)), cityname)
            gjson, colormap = choropleth.continuous_choropleth(filtered_data, factor=factor, layer_name = layer_name,
                                                                scale_name=scale_name, mid_color = mid_color, end_color = end_color,
                                                                show=False)
            colormap.add_to(basemap)
            gjson.add_to(basemap)
            BindColormap(gjson, colormap).add_to(basemap)

    # Add SF/MF permits
    sfpermits = pd.concat([austin_sf[['geometry']], dallas_sf[['geometry']], houston_sf[['geometry']]], axis = 0)
    sf_feature_group = FeatureGroup('SF Construction Permits, All 3 Cities', show = True)
    choropleth.make_marker_cluster(sfpermits, make_centroids = False, fast = True).add_to(sf_feature_group)
    sf_feature_group.add_to(basemap)

    mfpermits = pd.concat([austin_mf[['geometry']], dallas_mf[['geometry']], houston_mf[['geometry']]], axis = 0)
    mf_feature_group = FeatureGroup('MF Construction Permits, All 3 Cities', show = False)
    choropleth.make_marker_cluster(mfpermits, make_centroids = False, fast = True).add_to(mf_feature_group)
    mf_feature_group.add_to(basemap)


    # Put pricing data on the map, start with SF pricing data
    gjson, colormap = choropleth.continuous_choropleth(all_zip_geodata, factor = 'sfprice', layer_name = 'SF Prices in All Three Cities',
                                                       scale_name = 'Median SF Home Price ($)', show = False)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    # MF pricing data
    gjson, colormap = choropleth.continuous_choropleth(all_zip_geodata, factor = 'mfprice', layer_name = 'MF Prices in All Three Cities',
                                                       scale_name = 'Median MF Price (%)', mid_color = 'red', end_color = '#660000', show = False)
    colormap.add_to(basemap)
    gjson.add_to(basemap)
    BindColormap(gjson, colormap).add_to(basemap)

    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    folium.LayerControl().add_to(basemap)
    basemap.save('Figures/Bucket 2/PermitOutliers.html')

if income_histogram:

    calculate = False

    if calculate:

        # Get block data with income brackets
        block_data = sf.get_block_geodata(['X19_INCOME'], cities = ['Austin', 'Dallas'])
        austin_block_data = block_data['Austin']
        dallas_block_data = block_data['Dallas']
        factor_dictionary =   {'B19001e2':0, # Start values, the next value is the end of the bracket
                               'B19001e3':10000,
                               'B19001e4':15000,
                               'B19001e5':20000,
                               'B19001e6':25000,
                               'B19001e7':30000,
                               'B19001e8':35000,
                               'B19001e9':40000,
                               'B19001e10':45000,
                               'B19001e11':50000,
                               'B19001e12':60000,
                               'B19001e13':75000,
                               'B19001e14':100000,
                               'B19001e15':125000,
                               'B19001e16':150000,
                               'B19001e17':200000}

        data_features = [factor_dictionary[key] for key in factor_dictionary]
        austin_block_data = austin_block_data.rename(columns = factor_dictionary)
        dallas_block_data = dallas_block_data.rename(columns = factor_dictionary)

        # Now calculate averages by area
        dallas_zones = dallas_zones.loc[dallas_zones['broad_zone'].isin(['Single Family', 'Multifamily'])].to_crs({'init':'epsg:4326'})
        dallas_zones = sf.get_all_averages_by_area(dallas_block_data, dallas_zones, features = data_features,
                                 fillna = None)

        austin_zones = austin_zones.loc[austin_zones['broad_zone'].isin(['Single Family', 'Multifamily'])]
        austin_zones = sf.get_all_averages_by_area(austin_block_data, austin_zones, features = data_features,
                                 fillna = None)

        # Plot
        austin_zones['City'] = 'Austin'
        dallas_zones['City'] = 'Dallas'
        selected_columns = data_features
        selected_columns.extend(['broad_zone', 'City']) # We don't need geometry anymore
        all_zones = pd.concat([austin_zones[selected_columns], dallas_zones[selected_columns]], axis = 0)
        final_data = all_zones.groupby(['City', 'broad_zone']).sum()
        final_data.to_csv('Cached_Income_Data.csv')

    # Read in data just to get it in a consistent format - this is very fast anyways, it's like a 1 kb file
    final_data = pd.read_csv('Cached_Income_Data.csv')
    final_data = final_data.melt(var_name = 'Household_Income', value_name = 'Count', id_vars = ['City', 'broad_zone'])


    # Create categorical datatype for ordering
    from pandas.api.types import CategoricalDtype
    income_cat = CategoricalDtype(categories = final_data['Household_Income'].unique(), ordered = True)
    final_data['Household Income'] = final_data['Household_Income'].astype(income_cat)

    # Normalize by the total
    conditional_sums = final_data.groupby(['Household Income', 'broad_zone', 'City']).agg({'Count': 'sum'})
    final_data = conditional_sums.groupby(level = [1,2]).apply(lambda x: 100*x / x.sum()).reset_index()
    incomeplot = (ggplot(final_data, aes(x = 'Household Income', y = "Count", group = 'broad_zone', fill = 'broad_zone'))
               + geom_bar(stat="identity", position=position_dodge())
               + facet_wrap('~ City')
               + labs(title = 'Income by Base Residential Zone, Austin and Dallas',
                y = 'Percent of Households living in Housing Type')
               + theme(axis_text_x = element_text(rotation = 20, size = 8)))
    incomeplot.save(filename='Figures/Bucket 2/income_housing_typology.svg', width=15, height=8, bbox_inches='tight')


if construction_heatmap:

    # Get permit data ----------------------------------------

    # Get permit data - Dallas
    print('Processing Dallas permit data')
    dallas_permits = get_corrected_dallas_permit_data(path = dpm_save_path)
    dallas_permits = sf.process_points(dallas_permits)
    dallas_sf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Single')]
    dallas_mf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Multi')]

    # Get permit data - Austin
    austin_permits = process_austin_permit_data(searchfor=['101 single family houses',
                                                           '103 two family bldgs',
                                                           '104 three & four family bldgs',
                                                           '105 five or more family bldgs'],
                                                earliest=2013,
                                                permittypedesc='Building Permit',
                                                workclass='New')
    austin_permits = sf.process_points(austin_permits)
    austin_sf = austin_permits.loc[austin_permits['PermitClass'].str.contains('101 single')]
    austin_mf = austin_permits.loc[[not bool for bool in austin_permits['PermitClass'].str.contains('101 single')]]

    # Get permit data - Houston
    # Construction permit data
    print('Processing Houston permit data')
    houston_permit_data = process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE',
                                                                   'NEW AP', 'NEW HI-'],
                                                      searchin = ['PROJ_DESC'],
                                                      earliest = 2013, latest = None)

    # Subset to only include approved permits and nonempty geometries
    houston_permit_data = houston_permit_data.loc[houston_permit_data['Approval'] == 1.0]
    houston_permit_data = sf.process_points(houston_permit_data)
    houston_sf = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW S.F.',
                                                                                                 'NEW SF',
                                                                                                 'NEW TOWNHOUSE',
                                                                                                 'NEW SINGLE']))]
    houston_mf = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW AP',
                                                                                                 'NEW HI-']))]

    # Create heatmaps
    for city_name, city_input, sf_permits, mf_permits in zip(['Austin', 'Dallas', 'Houston'],
                                                             [austin_inputs, dallas_inputs, houston_inputs],
                                                             [austin_sf, dallas_sf, houston_sf],
                                                             [austin_mf, dallas_mf, houston_mf]):

        basemap = folium.Map([city_input.lat, city_input.long], zoom_start = 10)
        choropleth.heatmap(sf_permits,
                           name = 'Single Family Construction',
                           show = False,
                           radius = 13,
                           min_opacity = 0.5,
                           max_val = 1).add_to(basemap)
        choropleth.heatmap(mf_permits, name = 'Multifamily Construction',
                           show = False,
                           radius = 13,
                           min_opacity = 0.5).add_to(basemap)

        folium.TileLayer('cartodbdark_matter').add_to(basemap)
        folium.LayerControl().add_to(basemap)
        basemap.save('Figures/Bucket 2/Heatmaps/{}PermitHeatMap.html'.format(city_name))

