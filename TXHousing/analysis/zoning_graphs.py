"""Graphs which mostly rely on zoning data"""
import time
import shapely
import pandas as pd
import geopandas as gpd
from geopandas import plotting
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


def plot_broad_zones_dallas(colordic = {'Single Family':'#ff81c0', # Pink np.array((255, 129, 192), dtype = int)
                                        'Other Residential':'#c79fef', # Lavender np.array((199, 159, 239)), dtype = int)
                                        'Multifamily':'#840000', # Dark red np.array((132, 0, 0), dtype = int)
                                        'Other':'#96f97b'}, # Light green np.array((150, 249, 123), dtype = int)
                            save_path = 'Figures/Zoning/north_texas_base_zones.png',
                            minlat = 32.49, maxlat = 33.51, minlong = -97.71, maxlong = -96.29):
    """Plots an actual map of the closest broad_zones to the Dallas core. Colordic should map broad_zones to colors."""
    time0 = time.time()

    # Fig, ax
    fig, ax = plt.subplots()

    # Get subset of north texas data
    north_texas_data = zoning.north_texas_inputs.process_zoning_shapefile()
    north_texas_data['geometry'] = north_texas_data['geometry'].simplify(tolerance = 0.001)
    spatial_index = north_texas_data.sindex

    print('Subsetting')
    bounding_polygon = shapely.geometry.box(minlong, minlat, maxlong, maxlat)
    ids = list(spatial_index.intersection(bounding_polygon.bounds))
    subset = north_texas_data.iloc[ids]

    # Plot
    if colordic is not None:

        # Loop through zones with specific colors
        zones = subset['broad_zone'].unique()
        legend_handlers = []

        # Add zones
        for zone in zones:
            filtered_subset = subset.loc[subset['broad_zone'] == zone]
            filtered_subset.plot(ax = ax, color = colordic[zone], alpha = 0.6, label = zone)
            #gpd.plotting.plot_polygon_collection(ax, filtered_polygons, color = colordic[zone], alpha = 0.6, label = zone)
            legend_handlers.append(plt.scatter([], [], color =  colordic[zone]))

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title('Base Zoning Around Dallas, TX')
        ax.legend(tuple(legend_handlers), tuple(zones), fontsize = 6)

    else:
        subset.plot(ax = ax, column = 'broad_zone', legend = True, legend_kwds = {'fontsize':6}, cmap = 'Set1')

    # Plot and label counties
    counties = gpd.read_file(boundaries.county_boundaries_path)
    counties = counties.loc[(counties['NAME'].isin(['Dallas', 'Denton', 'Tarrant', 'Collin'])) & (counties['STATEFP'] == '48')]
    counties['coords'] = counties['geometry'].apply(lambda x: x.representative_point().coords[:][0])
    counties['geometry'] = counties['geometry'].apply(lambda x: x.boundary)
    counties.plot(ax = ax, facecolor = 'none', alpha = 0.5, edgecolor = 'black')
    for idx, row in counties.iterrows():
        ax.annotate(s=row['NAME'], xy=row['coords'], horizontalalignment='center')


    print('Saving')
    plt.savefig(save_path, dpi = 1000)
    print('Finished, took {}'.format(time.time() - time0))
