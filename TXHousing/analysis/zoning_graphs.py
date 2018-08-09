"""Graphs which mostly rely on zoning data"""
import time
import shapely
import pandas as pd
import geopandas as gpd
from .. import utilities
from ..data_processing import zoning, boundaries
from plotnine import *
import matplotlib.pyplot as plt

# Minimum lot size graph -------------------------------------
def plot_minimum_lot_size(savepath = 'Figures/Bucket 2/Minimum_Lot_Size_Residential_No_Agriculture.svg',
                          width = 10, height = 8):
    """Graphs average minimum lot size in Austin/Dallas.

    Methodological note: this treats zones with no minimum lot size as having a minimum lot size of zero. It also
    includes areas zoned as non-single family, which we might want to exclude in the future."""

    # Get zoning data and subset to exclude agricultural zones and nonresidential zones
    austin_zones = zoning.austin_inputs.process_zoning_shapefile(regulation_features = ['min_lot'])
    austin_zones = austin_zones.loc[(~austin_zones['base_zone'].str.contains('RR')) &
                                    (austin_zones['broad_zone'] != 'Other')]

    dallas_zones = zoning.dallas_inputs.process_zoning_shapefile(regulation_features = ['inv_density'])
    dallas_zones = dallas_zones.loc[(~dallas_zones['base_zone'].str.contains('A(A)', regex = False)) &
                                    (austin_zones['broad_zone'] != 'Other')]

    # Update min_lot/inv_density (which are the same thing) features to fill NaN values with zeros.
    austin_zones['min_lot'] = austin_zones['min_lot'].fillna(0)
    dallas_zones['inv_density'] = dallas_zones['inv_density'].fillna(0)

    # Polygons intersect rings calls

    austin_minlot_radii = utilities.measurements.polygons_intersect_rings(austin_zones, factor = 'min_lot',
                                                                          lat = zoning.austin_inputs.lat,
                                                                          long = zoning.austin_inputs.long,
                                                                          categorical = False,
                                                                          newproj = 'epsg:2277',
                                                                          step = 1,
                                                                          maximum = 15)

    dallas_minlot_radii = utilities.measurements.polygons_intersect_rings(dallas_zones, factor = 'inv_density',
                                                                          lat = zoning.dallas_inputs.lat,
                                                                          long = zoning.dallas_inputs.long,
                                                                          categorical = False,
                                                                          newproj = 'epsg:2276',
                                                                          step = 1,
                                                                          maximum = 15)


    minlotdata = pd.concat([dallas_minlot_radii, austin_minlot_radii], axis = 1)
    minlotdata.columns = ['Dallas', 'Austin']
    minlotdata['dist_to_center'] = minlotdata.index
    minlotdata = minlotdata.melt(id_vars = 'dist_to_center', value_name = 'min_lot', var_name = 'City')

    minlotplot = (ggplot(minlotdata, aes(x = 'dist_to_center', y = 'min_lot', fill = 'City'))
                  + geom_col(position = 'dodge', width = 0.9)
                  + theme_bw()
                  + scale_fill_manual(values=['#008000', 'cornflowerblue'])
                  + labs(title = 'Minimum Lot Size in Austin and Dallas, Excluding Agricultural and Nonresidential Land',
                         x = 'Distance from the City Center (Miles',
                         y = 'Average Lot Size (Square Feet)'))
    minlotplot.save(savepath, width = width, height = height)

def plot_hd_locations(save_path = 'Figures/Bucket 2/HDLocations.svg', width = 8, height = 5):
    """Graphs locations of historic districts in Austin, Dallas, Houston"""

    tx_hd_data = gpd.read_file(zoning.tx_hd_path).to_crs({'init':'epsg:4326'})

    # Get austin hds
    signature = '-HD'
    austin_nat_hds = tx_hd_data.loc[tx_hd_data['CITY'] == 'Austin']
    austin_zones = zoning.austin_inputs.process_zoning_shapefile()
    austin_local_hds =  austin_zones.loc[austin_zones[zoning.austin_inputs.feature].str.contains(signature)]

    # Get Dallas hds
    dallas_zones = zoning.dallas_inputs.process_zoning_shapefile()
    dallas_cds = dallas_zones.loc[dallas_zones['LONG_ZONE_'].apply(lambda x: x[0:2]) == 'CD']
    dallas_nat_hds = tx_hd_data.loc[tx_hd_data['CITY'] == 'Dallas']

    # Houston hds and others
    houston_local_hds = gpd.read_file(zoning.houston_historic_districts_path).to_crs({'init':'epsg:4326'})
    houston_nat_hds = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston']

    step = 0.5
    maximum = 5

    # Many polygons_intersect_rings calls
    austin_local_data = utilities.measurements.polygons_intersect_rings(austin_local_hds, factor=None,
                                                                        lat = zoning.austin_inputs.lat,
                                                                        long = zoning.austin_inputs.long,
                                                                        newproj='epsg:2277', step=step,
                                                                        maximum=maximum, group_outliers = False)
    austin_local_data.name = 'Austin Local Historic Districts'

    austin_nat_data = utilities.measurements.polygons_intersect_rings(austin_nat_hds, factor=None,
                                                                        lat = zoning.austin_inputs.lat,
                                                                        long = zoning.austin_inputs.long,
                                                                        newproj='epsg:2277', step=step,
                                                                        maximum=maximum, group_outliers = False)
    austin_nat_data.name = 'Austin National Historic Districts'

    dallas_cds_data = utilities.measurements.polygons_intersect_rings(dallas_cds, factor=None,
                                                                        lat = zoning.dallas_inputs.lat,
                                                                        long = zoning.dallas_inputs.long,
                                                                        newproj='epsg:2276', step=step,
                                                                        maximum=maximum, group_outliers = False)
    dallas_cds_data.name = 'Dallas Conservation Districts'
    dallas_nat_data = utilities.measurements.polygons_intersect_rings(dallas_nat_hds, factor=None,
                                                                    lat = zoning.dallas_inputs.lat,
                                                                    long = zoning.dallas_inputs.long,
                                                                    newproj='epsg:2276', step=step,
                                                                    maximum=maximum, group_outliers = False)
    dallas_nat_data.name = 'Dallas National Historic Districts'

    houston_local_data = utilities.measurements.polygons_intersect_rings(houston_local_hds, factor=None,
                                                                        lat = zoning.houston_inputs.lat,
                                                                        long = zoning.houston_inputs.long,
                                                                        newproj='epsg:2278', step=step,
                                                                        maximum=maximum, group_outliers = False)
    houston_local_data.name = 'Houston Local Historic Districts'


    houston_nat_data = utilities.measurements.polygons_intersect_rings(houston_nat_hds, factor=None,
                                                                        lat = zoning.houston_inputs.lat,
                                                                        long = zoning.houston_inputs.long,
                                                                        newproj='epsg:2278', step=step,
                                                                        maximum=maximum, group_outliers = False)
    houston_nat_data.name = 'Houston National Historic Districts'

    # Combine and plot
    all_loc_data = 100 * pd.concat(
        [austin_local_data, austin_nat_data, dallas_cds_data, dallas_nat_data, houston_local_data, houston_nat_data],
        axis=1)
    all_loc_data['dist'] = all_loc_data.index
    all_loc_data = pd.melt(all_loc_data, var_name='type', value_name='percent', id_vars=['dist'])
    all_loc_data['city'] = all_loc_data['type'].apply(lambda x: x.split(' ')[0])
    all_loc_data['District Type'] = all_loc_data['type'].apply(lambda x: x.split(' ')[1])

    histlocations = (ggplot(all_loc_data, aes(x='dist', y='percent', group='District Type', fill='District Type'))
                     + geom_col(position='dodge')
                     + facet_wrap('~ city')
                     + labs(title = 'Locations of Historic Districts in Austin, Dallas, and Houston',
                            x = 'Distance from City Center (Miles)',
                            y = 'Percent of Land in Historic Districts')
                     + theme_bw())
    histlocations.save(save_path, width=width, height=height)

def plot_broad_zones_proportion():
    """Plot proportion of broad_zones by distance from city center, excluding nonresidential and agricultural land."""

    # Get zoning data and subset to exclude agricultural zones and nonresidential zones
    austin_zones = zoning.austin_inputs.process_zoning_shapefile(regulation_features = ['min_lot'])
    austin_zones = austin_zones.loc[(~austin_zones['base_zone'].str.contains('RR')) &
                                    (austin_zones['broad_zone'] != 'Other')]

    dallas_zones = zoning.dallas_inputs.process_zoning_shapefile(regulation_features = ['inv_density'])
    dallas_zones = dallas_zones.loc[(~dallas_zones['base_zone'].str.contains('A(A)', regex = False)) &
                                    (dallas_zones['broad_zone'] != 'Other')]

    dallas_zone_rings = utilities.measurements.polygons_intersect_rings(dallas_zones, factor = 'broad_zone',
                                                                        lat = zoning.dallas_inputs.lat,
                                                                        long = zoning.dallas_inputs.long,
                                                                        categorical = True,
                                                                        newproj='epsg:2276', step=1, maximum = 10)
    dallas_zone_rings['City'] = "Dallas"

    austin_zone_rings = utilities.measurements.polygons_intersect_rings(austin_zones, factor = 'broad_zone',
                                                                        lat = zoning.austin_inputs.lat,
                                                                        long = zoning.austin_inputs.long,
                                                                        categorical = True,
                                                                        newproj='epsg:2277', step=1, maximum = 10)
    austin_zone_rings['City'] = 'Austin'

    # Combine and divide by totals
    zone_rings = pd.concat([dallas_zone_rings, austin_zone_rings], axis = 0)
    zone_rings['dist_to_center'] = zone_rings.index
    zone_rings = zone_rings.melt(id_vars = ['City', 'dist_to_center'], var_name = 'broad_zone', value_name = 'percent')

    # Plot Single Family
    sfrings = zone_rings.loc[zone_rings['broad_zone'] == 'Single Family']
    sfringsplot = (ggplot(sfrings, aes(x = 'dist_to_center', y = 'percent', fill = 'City'))
                   + geom_col(position = 'dodge', width = 0.6)
                   + theme_bw()
                   + scale_fill_manual(['cornflowerblue', '#8FBC8F'])
                   + labs(title = 'Single Family Zoning in Austin and Dallas, Excluding Agricultural and Nonresidential Land',
                          x = 'Distance from the city center (Miles)',
                          y = 'Percentage of Residential Land Zoned'))
    sfringsplot.save('Figures/Bucket 2/SFZoningRings.svg', width = 8, height = 5)

    # Plot Multifamily
    mfrings = zone_rings.loc[zone_rings['broad_zone'] == 'Multifamily']
    mfringsplot = (ggplot(mfrings, aes(x = 'dist_to_center', y = 'percent', fill = 'City'))
                   + geom_col(position = 'dodge', width = 0.6)
                   + theme_bw()
                   + scale_fill_manual(['cornflowerblue', '#8FBC8F'])
                   + labs(title = 'Multifamily Zoning in Austin and Dallas, Excluding Agricultural and Nonresidential Land',
                          x = 'Distance from the city center (Miles)',
                          y = 'Percentage of Residential Land Zoned'))
    mfringsplot.save('Figures/Bucket 2/MFZoningRings.svg', width = 8, height = 5)

# Matplotlib helper functions
def add_counties(ax, county_list):
    """ Given a matplotlib axis, adds texas county outlines to it"""

    # Plot and label counties
    counties = gpd.read_file(boundaries.county_boundaries_path)
    counties = counties.loc[(counties['NAME'].isin(county_list)) & (counties['STATEFP'] == '48')]
    counties['coords'] = counties['geometry'].apply(lambda x: x.representative_point().coords[:][0])
    counties['geometry'] = counties['geometry'].apply(lambda x: x.boundary)
    counties.plot(ax = ax, facecolor = 'none', alpha = 0.5, edgecolor = 'black')
    for idx, row in counties.iterrows():
        ax.annotate(s=row['NAME'], xy=row['coords'], horizontalalignment='center')

def map_broad_zones(data, county_list, name, minlat, maxlat, minlong, maxlong, save_path, colordic = None):
    """Plots an actual map of broad_zones in a data source within the lat/long bounds. If you want to do this for Austin
    or Dallas, just use the map_broad_zones_dallas_austin wrapper.

    :param data: A gdf including a geometry and base_zone column.
    :param county_list: A list of county outlines to plot on the data.
    :param name: A name to be used in titling the egraph.
    :param minlat, maxlat, minlong, maxlong: Bounds of the graph
    :param save_path: A path at which to save the map.
    :param colordic: Dictionary which maps broad_zones to colors."""

    time0 = time.time()

    # Fig, ax
    fig, ax = plt.subplots()

    # Get subset of north texas data
    data['geometry'] = data['geometry'].simplify(tolerance = 0.001)
    spatial_index = data.sindex

    print('Subsetting')
    bounding_polygon = shapely.geometry.box(minlong, minlat, maxlong, maxlat)
    ids = list(spatial_index.intersection(bounding_polygon.bounds))
    subset = data.iloc[ids]

    # Plot
    if colordic != 'none':

        if colordic is None:
            colordic = {'Single Family': '#ff81c0',  # Pink np.array((255, 129, 192), dtype = int)
                         'Other Residential': '#c79fef',  # Lavender np.array((199, 159, 239)), dtype = int)
                         'Multifamily': '#840000',  # Dark red np.array((132, 0, 0), dtype = int)
                         'Other': '#96f97b'}

        # Loop through zones with specific colors
        zones = subset['broad_zone'].unique()
        legend_handlers = []

        # Add zones
        for zone in zones:
            filtered_subset = subset.loc[subset['broad_zone'] == zone]
            filtered_subset.plot(ax = ax, color = colordic[zone], alpha = 0.6, label = zone)
            legend_handlers.append(plt.scatter([], [], color =  colordic[zone]))

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_ylim(minlat, maxlat)
        ax.set_xlim(minlong, maxlong)
        ax.set_title('Base Zoning Around {}, TX'.format(name))
        ax.legend(tuple(legend_handlers), tuple(zones), fontsize = 6)

    else:
        subset.plot(ax = ax, column = 'broad_zone', legend = True, legend_kwds = {'fontsize':6}, cmap = 'Set1')

    add_counties(ax, county_list = county_list)

    print('Saving')
    plt.savefig(save_path, dpi = 1000)
    print('Finished, took {}'.format(time.time() - time0))

def map_broad_zones_dallas_austin(plot_austin = True, plot_dallas = True):
    """ Wrapper for map_broad_zones, just plots them around austin/dallas. """

    if plot_austin:
        austin_data = zoning.get_austin_surrounding_zones()
        map_broad_zones(data = austin_data,
                        county_list = ['Travis', 'Williamson'],
                        minlat=30.09, maxlat=30.76, minlong=-98.11, maxlong=-97.43,
                        name = 'Austin',
                        save_path = 'Figures/Zoning/Austin_base_zones.png')

    if plot_dallas:
        north_texas_data = zoning.north_texas_inputs.process_zoning_shapefile()
        map_broad_zones(data = north_texas_data,
                        county_list = ['Dallas', 'Denton', 'Tarrant', 'Collin'],
                        minlat=32.49, maxlat=33.51, minlong=-97.71, maxlong=-96.29,
                        name = 'Dallas',
                        save_path = 'Figures/Zoning/north_texas_base_zones.png')



def plot_zone_income_histogram(save_path = 'Figures/Zoning/income_housing_typology.svg',
                               cache_path = 'shared_data/calculations/zoning_income_data.csv'):
    """ Plots the distribution of incomes conditional on broad zones in Austin and Dallas"""

    # Get block data with income brackets
    block_data = boundaries.BlockBoundaries(['X19_INCOME'], cities = ['Austin', 'Dallas'])
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
    block_data.data = block_data.data.rename(columns = factor_dictionary)

    # Austin/Dallas Zones
    dallas_zones = zoning.dallas_inputs.process_zoning_shapefile()
    dallas_zones = dallas_zones.loc[dallas_zones['broad_zone'].isin(['Single Family', 'Multifamily'])]
    austin_zones = zoning.austin_inputs.process_zoning_shapefile()
    austin_zones = austin_zones.loc[austin_zones['broad_zone'].isin(['Single Family', 'Multifamily'])]

    # Pull income statistics
    austin_zones = block_data.push_features(austin_zones, data_features)
    dallas_zones = block_data.push_features(dallas_zones, data_features)

    # Plot
    austin_zones['City'] = 'Austin'
    dallas_zones['City'] = 'Dallas'
    selected_columns = data_features
    selected_columns.extend(['broad_zone', 'City']) # We don't need geometry anymore
    all_zones = pd.concat([austin_zones[selected_columns], dallas_zones[selected_columns]], axis = 0)
    final_data = all_zones.groupby(['City', 'broad_zone']).sum()
    final_data.to_csv(cache_path)

    # Read in data just to get it in a consistent format - this is very fast anyways, it's like a 1 kb file
    final_data = pd.read_csv(cache_path)
    final_data = final_data.melt(var_name = 'Household_Income', value_name = 'Count', id_vars = ['City', 'broad_zone'])

    final_data = utilities.measurements.order_radii(final_data, feature = 'Household_Income')

    # Normalize by the total
    conditional_sums = final_data.groupby(['Household_Income', 'broad_zone', 'City']).agg({'Count': 'sum'})
    final_data = conditional_sums.groupby(level = [1,2]).apply(lambda x: 100*x / x.sum()).reset_index()

    incomeplot = (ggplot(final_data, aes(x = 'Household_Income', y = "Count", group = 'broad_zone', fill = 'broad_zone'))
               + geom_bar(stat="identity", position=position_dodge())
               + facet_wrap('~ City')
               + theme_bw()
               + labs(title = 'Income by Base Residential Zone, Austin and Dallas',
                      x = 'Household Income Bracket',
                      y = 'Percent of SF/MF Households in Income Bracket Given Income')
               + theme(axis_text_x = element_text(rotation = 20, size = 8)))
    incomeplot.save(filename=save_path, width=15, height=8, bbox_inches='tight')
