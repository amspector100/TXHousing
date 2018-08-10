"""Graphs which mostly rely on permit data"""

import time
import pandas as pd
from .. import utilities
from ..data_processing import permit, boundaries, property, zoning
from plotnine import *

def plot_permit_scatterplot():
    """Scatterplots of the number of construction of permits in a variety of zipcodes against the median housing price
     of the zipcode."""

    time0 = time.time()

    # Get permit data - Dallas
    print('Getting Dallas permit data, time is {}'.format(time.time() - time0))
    dallas_permits = permit.get_corrected_dallas_permit_data()
    dallas_permits = utilities.simple.process_points(dallas_permits)
    dallas_sf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Single')]
    dallas_mf = dallas_permits.loc[dallas_permits['Permit Type'].str.contains('Multi')]

    # Get permit data - Austin
    print('Getting Austin permit data, time is {}'.format(time.time() - time0))
    austin_permits = permit.process_austin_permit_data(searchfor=['101 single family houses',
                                                           '103 two family bldgs',
                                                           '104 three & four family bldgs',
                                                           '105 five or more family bldgs'],
                                                earliest=2013,
                                                permittypedesc='Building Permit',
                                                workclass='New')
    austin_permits = utilities.simple.process_points(austin_permits)
    austin_sf = austin_permits.loc[austin_permits['PermitClass'].str.contains('101')]
    austin_mf = austin_permits.loc[~austin_permits['PermitClass'].str.contains('101')]

    # Get permit (construction) data - Houston
    print('Getting Houston permit data, time is {}'.format(time.time() - time0))
    houston_permit_data = permit.process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE',
                                                                          'NEW TOWNHOUSE', 'NEW AP', 'NEW HI-'],
                                                              searchin = ['PROJ_DESC'],
                                                              earliest = 2013, latest = None)

    # Subset to only include approved permits and nonempty geometries
    houston_permit_data = houston_permit_data.loc[houston_permit_data['Approval'] == 1.0]
    houston_permit_data = utilities.simple.process_points(houston_permit_data)
    houston_sf = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW S.F.',
                                                                                                 'NEW SF',
                                                                                                 'NEW TOWNHOUSE',
                                                                                                 'NEW SINGLE']))]
    houston_mf = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW AP',
                                                                                                 'NEW HI-']))]

    # Now get zip data and merge with realtor data
    print('Working with Zip Boundaries and Permit Data, time is {}'.format(time.time() - time0))
    zipdata = boundaries.ZipBoundaries()
    zipdata.add_property_data(property.realtor_core_inventory_sf, features = ['Median Listing Price'], rsuffix = '_sf')
    zipdata.add_property_data(property.realtor_core_inventory_mf, features = ['Median Listing Price'], rsuffix = '_mf')
    zipdata.data['City'] = float("NaN")

    for ziplist, cityname in zip([boundaries.austin_zips, boundaries.dallas_zips, boundaries.houston_zips],
                                 ['Austin', 'Dallas', 'Houston']):
        ziplist = [z for z in ziplist if z in zipdata.data.index]
        zipdata.data.loc[ziplist, 'City'] = cityname

    zipdata.data = zipdata.data.loc[zipdata.data['City'].notnull()]

    # Spatial join with austin/dallas permit data - houston permit data already has zipcode information in it
    for pdata, colname in zip([austin_sf, austin_mf, dallas_sf, dallas_mf],
                              ['austin_sf', 'austin_mf', 'dallas_sf', 'dallas_mf']):
        intersections = zipdata.fast_intersection(pdata)
        counted = intersections.groupby(intersections).count()
        zipdata.data[colname] = counted

    # Houston sf/mf
    zipdata.data['houston_sf'] = houston_sf.groupby(['Zipcode'])['geometry'].count()
    zipdata.data['houston_mf'] = houston_mf.groupby(['Zipcode'])['geometry'].count()

    # Start to melt and get ready to plot for single family
    sfresult = zipdata.data[['austin_sf', 'dallas_sf', 'houston_sf', 'Median Listing Price', 'City']]
    sfresult = sfresult.melt(id_vars = ['City', 'Median Listing Price'], var_name = 'permit_type', value_name = 'num_permits')

    y_max = 1000
    sfplot = (ggplot(sfresult, aes(x = 'Median Listing Price', y = 'num_permits', color = 'City'))
              + geom_point()
              + stat_smooth(method='lowess', span=0.5)
              + facet_wrap('City', scales='fixed')
              + scale_x_continuous(limits=(0, 2200000), breaks=(0, 1000000, 2000000))
              + scale_y_continuous(limits=(0, y_max))
              + labs(title = 'New Single Family Construction Permits against SF Housing Price',
                     caption = 'Permit Data from the Cities of Austin, Dallas, and Houston /n, Pricing data from Realtor')
              + xlab('Median Single Family Home Price (by Zip Code)')
              + ylab('Number of SF Construction Permits Issued In Past 5 Years'))

    # Start to melt and get ready to plot for multifamily
    mfresult = zipdata.data[['austin_mf', 'dallas_mf', 'houston_mf', 'Median Listing Price_mf', 'City']]
    mfresult = mfresult.melt(id_vars = ['City', 'Median Listing Price_mf'], var_name = 'permit_type', value_name = 'num_permits').drop('permit_type', axis = 1)

    mfplot = (ggplot(mfresult, aes(x = 'Median Listing Price_mf', y = 'num_permits', color = 'City'))
          + geom_point()
          + stat_smooth(method = 'lowess', span = 0.5)
          + facet_wrap('City', scales='fixed'))

    width = 10
    height = 5
    sfplot.save(filename='Figures/Permit/construction_scatter_sf.svg', width=width, height=height, bbox_inches = 'tight')
    mfplot.save(filename='Figures/Permit/construction_scatter_mf.svg', width=width, height=height, bbox_inches = 'tight')


def plot_permit_locations(save_path = 'Figures/Permit/all_cities_permit_rings.svg', width = 10, height = 8, step = 1,
                          maximum = 10):
    """ Plots the number of permits issued per square mile, conditional on from city center & broad_zone"""

    time0 = time.time()

    # Austin
    austin_permit_data = permit.process_austin_permit_data(searchfor=['101 single family houses',
                                                                      '103 two family bldgs',
                                                                      '104 three & four family bldgs',
                                                                      '105 five or more family bldgs'],
                                                           earliest=2013)
    austin_permit_data['broad_zone'] = austin_permit_data['PermitClass'].apply(lambda x: 'SF' if '101' in x else 'MF')
    austin_permit_rings = utilities.measurements.points_intersect_rings(austin_permit_data,
                                                                        lat = zoning.austin_inputs.lat,
                                                                        long = zoning.austin_inputs.long,
                                                                        factor = 'broad_zone', step = step,
                                                                        categorical = True, maximum = maximum)
    austin_permit_rings['City'] = 'Austin'
    print('Finished with Austin, time is {}'.format(time.time() - time0))

    # Dallas
    dallas_permit_data = permit.get_corrected_dallas_permit_data()  # Get dallas permit data from the saved file with CORRECTED geocodes and the original columns.
    dallas_permit_data['broad_zone'] = dallas_permit_data['Permit Type'].apply(lambda x: 'SF' if 'Single F' in x else 'MF')
    dallas_permit_rings = utilities.measurements.points_intersect_rings(dallas_permit_data,
                                                                        lat = zoning.dallas_inputs.lat,
                                                                        long = zoning.dallas_inputs.long,
                                                                        factor='broad_zone', step=step,
                                                                        categorical=True, maximum=maximum)
    dallas_permit_rings['City'] = 'Dallas'
    print('Finished with Dallas, time is {}'.format(time.time() - time0))

    houston_permit_data = permit.process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE',
                                                                          'NEW AP', 'NEW HI-'],
                                                             searchin = ['PROJ_DESC'],
                                                             earliest = 2013, latest = None)
    houston_permit_data['broad_zone'] = houston_permit_data['PROJ_DESC'].apply(lambda x: 'SF' if ('NEW S' in x or 'NEW T' in x) else 'MF')

    houston_permit_rings = utilities.measurements.points_intersect_rings(houston_permit_data,
                                                                         lat = zoning.houston_inputs.lat,
                                                                         long = zoning.houston_inputs.long,
                                                                         factor='broad_zone', step=step,
                                                                         categorical=True, maximum=maximum)
    houston_permit_rings['City'] = 'Houston'
    print('Finished with Houston, time is {}'.format(time.time() - time0))

    all_permit_data = pd.concat([austin_permit_rings, dallas_permit_rings, houston_permit_rings], axis = 0)
    all_permit_data['dist_to_center'] = all_permit_data.index

    all_permit_data = pd.melt(all_permit_data, var_name = 'permit_type', value_name = 'num_permits', id_vars = ['City', 'dist_to_center'])
    all_permit_data['permit_type'] = all_permit_data['permit_type'].apply(lambda x: 'Single Family' if x == 'SF' else 'Multifamily')

    permit_plot = (ggplot(all_permit_data, aes(x = 'dist_to_center', y = 'num_permits', fill = 'City'))
                   + geom_col(position = 'dodge')
                   + facet_wrap('~permit_type', scales = 'free')
                   + theme_bw()
                   + labs(title = 'Construction Permits Issued by Distance from City Center in Texas Triangle',
                          x = 'Distance from City Center (Miles)',
                          y = 'Number of Construction Permits Issued per Square Mile')
                   + scale_x_discrete(expand=(0, 0.1)))
    permit_plot.save(save_path, width = width, height = height)
