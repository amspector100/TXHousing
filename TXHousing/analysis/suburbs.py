""" This file contains most of the graphing for the suburbs analysis. It does not rely on raw parcel data, but the
cached versions instead."""

import numpy as np
import pandas as pd
import geopandas as gpd
import folium
from plotnine import *

from ..utilities import measurements
from ..data_processing import boundaries, zoning, parcel
from . import choropleth

# Functions which find paths - their names say it all
def get_municipality_choropleth_path(name):
    return 'Figures/Suburbs/{}_suburb_choropleth.html'.format(name)

def land_use_by_municipality_path(name):
    return 'shared_data/calculations/{}_land_use_by_municipality.csv'.format(name)

def lot_size_by_municipality_path(name):
    return 'shared_data/calculations/{}_lot_size_by_municipality.csv'.format(name)

# Plots percent of workers
def texas_job_centers(num_polygons = 4000):
    """ Choropleth of the percent of workers who work in their place of residence for Austin, Dallas, Houston """

    # Get block geodata for all of Texas and calculate features
    block_data = boundaries.BlockBoundaries(['X08_COMMUTING', 'X01_AGE_AND_SEX'], cities = None)

    # Female and male workers working in their place of residence
    block_data.data['local_workers'] = block_data.data['B08008e3'] + block_data.data['B08008e8']
    # Total number of male and female workers in the block group
    block_data.data['total_workers'] = block_data.data['B08008e2'] + block_data.data['B08008e7']
    block_data.data['local_workers_pct'] = 100*block_data.data['local_workers'].divide(block_data.data['total_workers']).fillna(0)

    # Get sindex and loop through blocks around austin/dallas/houston
    spatial_index = block_data.sindex

    for name, zoning_input in zip(['Austin', 'Dallas', 'Houston'],
                                  [zoning.austin_inputs, zoning.dallas_inputs, zoning.houston_inputs]):

        # Get basemap
        basemap = folium.Map([zoning_input.lat, zoning_input.long], zoom_start=10)

        # Query and find nearest neighbors, subset
        nearest_index = list(spatial_index.nearest((zoning_input.long, zoning_input.lat), num_results = num_polygons))
        city_data = block_data.data.iloc[nearest_index]

        # Graph
        choropleth.continuous_choropleth(city_data, factor = 'local_workers_pct',
                                         layer_name='Percent of Workers Working in Place of Residence',
                                         colors = ['white', 'green', 'blue'],
                                         quants = None,
                                         method = 'linear',
                                         show = False, basemap = basemap)
        # Save
        folium.TileLayer('cartodbdark_matter').add_to(basemap)
        folium.LayerControl().add_to(basemap)
        basemap.save('Figures/Suburbs/{}_job_choropleth.html'.format(name))
        print('Graphed for {}'.format(name))

    print('Finished')


def analyze_land_use_by_metro(name):
    """ Calculates mean and median lot size as well as the land use by municipality around Austin, Dallas, or Houston,
    then caches them in shared_data/calculations."""

    # Read data
    path = parcel.get_cached_all_parcel_path_csv(name)
    data = pd.read_csv(path, engine = 'python')

    # Calculate percent of area zoned
    zone_areas = data.groupby(['broad_zone', 'place'])['area_sqft'].sum()
    municipality_areas = data.groupby(['place'])['area_sqft'].sum()
    final_data = zone_areas.divide(municipality_areas)
    final_data = final_data.reset_index()
    final_data.columns = ['broad_zone', 'place', 'Percent of Land Used']
    final_data['Percent of Land Used'] = 100*final_data['Percent of Land Used']
    final_data.to_csv(land_use_by_municipality_path(name))

    # Calculate lot sizes (average and mean)
    counts = data.groupby(['broad_zone', 'place'])['area_sqft'].count()
    lotsize_means = zone_areas.divide(counts).reset_index()
    lotsize_means.columns = ['broad_zone', 'place', 'mean_lot_size']
    lotsize_medians = data.groupby(['broad_zone', 'place'])['area_sqft'].median().reset_index()
    lotsize_medians.columns = ['broad_zone', 'place', 'median_lot_size']
    lotsize_calculations = lotsize_means.merge(lotsize_medians, on = ['broad_zone', 'place'], how = 'outer')
    lotsize_calculations.to_csv(lot_size_by_municipality_path(name))

def suburbs_choropleth(name, zoning_input, colors = ['purple', 'red', 'orange', 'yellow', 'white'], method = 'linear'):
    """ Folium choropleths of land use, lotsize, and more by municipality in the suburbs. """

    # Read lotsize, landuse data
    lotsize_data = pd.read_csv(lot_size_by_municipality_path(name), index_col = 0)
    landuse_data = pd.read_csv(land_use_by_municipality_path(name), index_col = 0)

    # Join to polygon data. The dist_to_center calculations are for sorting.
    texas_place_shapes = gpd.read_file(boundaries.texas_places_path).to_crs({'init':'epsg:4326'})
    texas_place_shapes['dist_to_center'] = measurements.calculate_dist_to_center(texas_place_shapes,
                                                                                 lat = zoning_input.lat,
                                                                                 long = zoning_input.long)
    texas_place_shapes = texas_place_shapes.sort_values('dist_to_center')
    texas_place_shapes = texas_place_shapes.drop_duplicates(subset = ['NAME'], keep = 'first')
    landuse_data = texas_place_shapes.merge(landuse_data, left_on = 'NAME', right_on = 'place', how = 'right')
    lotsize_data = texas_place_shapes.merge(lotsize_data, left_on = 'NAME', right_on = 'place', how = 'right')

    # Plot
    basemap = folium.Map([zoning_input.lat, zoning_input.long], zoom_start = 9)

    # Lot size for single family only
    choropleth.continuous_choropleth(lotsize_data.loc[lotsize_data['broad_zone'] == 'Single Family'],
                                     'mean_lot_size',
                                     layer_name = 'Mean Single Family Lot Size',
                                     scale_name = 'Mean Single Family Lot Size (Sqft)',
                                     quants = [1/10, 2/10, 3/10, 5/10, 8/10],
                                     colors=colors, method=method, round_method = 'log10', show=False,
                                     basemap=basemap)

    # Percent of land used for single family
    choropleth.continuous_choropleth(landuse_data.loc[landuse_data['broad_zone'] == 'Single Family'],
                                     'Percent of Land Used',
                                     layer_name = 'Percent of Land Developed as Single Family',
                                     scale_name = 'Percent of Land Developed as Single Family',
                                     quants = None,
                                     colors=colors, method=method, round_method = 'int', show=False,
                                     basemap=basemap)

    # Percent of land used for Multifamily
    choropleth.continuous_choropleth(landuse_data.loc[landuse_data['broad_zone'] == 'Multifamily'],
                                     'Percent of Land Used',
                                     layer_name = 'Percent of Land Developed as Multifamily',
                                     scale_name = 'Percent of Land Developed as Multifamily',
                                     quants = None,
                                     colors=colors, method=method, round_method = 'int', show=False,
                                     basemap=basemap)

    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    folium.TileLayer('CartoDB positron').add_to(basemap)
    folium.LayerControl().add_to(basemap)
    basemap.save(get_municipality_choropleth_path(name))

# Globals for the suburbs scatterplot

dallas_urban_areas = ['Dallas--Fort Worth--Arlington, TX', 'Denton--Lewisville, TX', 'McKinney, TX']
dallas_job_centers = ['Plano', 'Irving', 'Fort Worth', 'Arlington', 'Lewisville', 'McKinney', 'Rockwall', 'Garland',
                      'Denton', 'Frisco']

houston_urban_areas = ['Houston, TX'] # 'Conroe--The Woodlands, TX', 'Texas City, TX' could be but are not included
houston_job_centers = ['Rosenberg', 'Sugarland', 'The Woodlands', 'Katy', 'Pearland', 'La Porte', 'Friendswood']

austin_urban_areas = ['Austin, TX']
austin_job_centers = ['Round Rock', 'Georgetown', 'Cedar Park', 'Leandor', 'Taylor', 'Elgin', 'Bastrop', 'Lakeway']

def suburbs_scatterplot():

    # Get place shapes
    texas_place_shapes = gpd.read_file(boundaries.texas_places_path)
    texas_place_shapes = measurements.get_area_in_units(texas_place_shapes, final_projection={'init':'epsg:4326'})

    # Subset to only include place shapes inside the predefined uas for each city
    ua_shapes = boundaries.Boundaries(boundaries.ua_path)
    ua_shapes.data = ua_shapes.data.loc[ua_shapes.data['NAME10'].str.contains('TX')]
    mapper = ua_shapes.fast_intersection(texas_place_shapes)
    texas_place_shapes['ua'] = texas_place_shapes.index.map(mapper).map(ua_shapes.data['NAME10'])
    all_uas = austin_urban_areas + dallas_urban_areas + houston_urban_areas
    texas_place_shapes = texas_place_shapes.loc[texas_place_shapes['ua'].isin(all_uas)]

    # Iterate through and get some geospatial features and data
    names = ['austin', 'dallas', 'houston']
    inputs = [zoning.austin_inputs, zoning.dallas_inputs, zoning.houston_inputs]
    jobcenter_lists = [austin_job_centers, dallas_job_centers, houston_job_centers]
    lotsize_data = pd.DataFrame()
    landuse_data = pd.DataFrame()

    for name, zoning_input, jobcenters in zip(names, inputs, jobcenter_lists):

        # Geodata
        new_shapes = texas_place_shapes.copy()
        new_shapes['dist_to_center'] = measurements.calculate_dist_to_center(new_shapes,
                                                                             lat = zoning_input.lat,
                                                                             long = zoning_input.long)
        new_shapes = new_shapes.sort_values('dist_to_center')
        new_shapes = new_shapes.drop_duplicates(subset = 'NAME', keep = 'first')
        new_shapes = new_shapes[['NAME', 'area', 'dist_to_center']]

        # Landuse by municipality
        new_landuse_data = pd.read_csv(land_use_by_municipality_path(name), index_col = 0)
        new_landuse_data['city'] = name
        new_landuse_data['jobcenter'] = new_landuse_data['place'].apply(lambda x: x in jobcenters)
        new_landuse_data = new_shapes.merge(new_landuse_data, left_on = 'NAME', right_on = 'place', how = 'right')
        landuse_data = landuse_data.append(new_landuse_data)

        # Lotsize by municipality
        new_lotsize_data = pd.read_csv(lot_size_by_municipality_path(name), index_col = 0)
        new_lotsize_data['city'] = name
        new_lotsize_data['jobcenter'] = new_lotsize_data['place'].apply(lambda x: x in jobcenters)
        new_lotsize_data = new_shapes.merge(new_lotsize_data, left_on = 'NAME', right_on = 'place', how = 'right')
        lotsize_data = lotsize_data.append(new_lotsize_data)

    # This data is by city
    landuse_data = landuse_data.loc[landuse_data['dist_to_center'] < 50]
    # Get rid of huge outliers
    lotsize_data = lotsize_data.loc[(lotsize_data['dist_to_center'] < 50) & (lotsize_data['mean_lot_size'] < 150000)] #

    sfplot = (ggplot(landuse_data.loc[landuse_data['broad_zone'] == 'Single Family'],
                    aes(x = 'dist_to_center', y = 'Percent of Land Used', size = 'area',
                        shape = 'jobcenter', fill = 'city', label = 'place'))
             + geom_point()
             + geom_text(nudge_y = 0.025) # 0.01 for unfacceted version
             + facet_wrap('~city')
             + labs(title = 'Percent of Land Developed as Single Family in Texas Triangle Urban Areas',
                    x = 'Distance from City Center (Miles)',
                    y = 'Percent of Land Used as Single Family'))
    sfplot.save('Figures/Suburbs/sf_landuse_scatterplot.svg', width = 8.5, height = 11)

    lplot = (ggplot(lotsize_data.loc[lotsize_data['broad_zone'] == 'Single Family'],
                    aes(x = 'dist_to_center', y = 'mean_lot_size', size = 'area', shape = 'jobcenter', fill = 'city', label = 'place'))
             + geom_point()
             + scale_y_log10()
             + geom_text(nudge_y=0.025)
             + facet_wrap('~city')
             + labs(title='Average Single Family Lotsize in Texas Triangle Urban Areas',
                    x='Distance from City Center (Miles)',
                    y='Average Single Family Lotsize'))
    lplot.save('Figures/Suburbs/sf_lotsize_scatterplot.svg', width = 8.5, height = 11)



# Commuting ------------------------------------------------------------------------------------------------------
def analyze_transportation_networks(names, zoning_inputs, num_blocks = 5000, step = 2.5, maximum = 60):

    block_data = boundaries.BlockBoundaries(['X08_COMMUTING', 'X01_AGE_AND_SEX'], cities = None)

    # Get total/local workers

    # Female and male workers working in their place of residence
    block_data.data['local_workers'] = block_data.data['B08008e3'] + block_data.data['B08008e8']
    # Total number of male and female workers in the block group
    block_data.data['total_workers'] = block_data.data['B08008e2'] + block_data.data['B08008e7']

    # Create spatial indexes and subset
    spatial_index = block_data.data.sindex
    all_data = gpd.GeoDataFrame()
    for name, zoning_input in zip(names, zoning_inputs):

        # Query and find nearest neighbors, subset
        nearest_index = list(spatial_index.nearest((zoning_input.long, zoning_input.lat), num_results=num_blocks))
        city_data = block_data.data.iloc[nearest_index]
        city_data['dist_to_center'] = measurements.calculate_dist_to_center(city_data,
                                                                            lat = zoning_input.lat,
                                                                            long = zoning_input.long)
        city_data['City'] = name
        all_data = pd.concat([all_data, city_data])

    # Average commute time

    all_data['smoothed_dist_to_center'] =  all_data['dist_to_center'].apply(lambda x: step*(np.round(x/step)))
    commute_brackets = {'B08303e2':2.5,
                        'B08303e3':7.5,
                        'B08303e4':12.5,
                        'B08303e5':17.5,
                        'B08303e6':22.5,
                        'B08303e7':27.5,
                        'B08303e8':32.5,
                        'B08303e9':37.5,
                        'B08303e10':42.5,
                        'B08303e11':52.5,
                        'B08303e12':75,
                        'B08303e13':90}
    all_data['total_commute_time'] = 0
    for key in commute_brackets:
        all_data['total_commute_time'] += all_data[key]*commute_brackets[key]
    print(all_data[['total_commute_time', 'B08135e1']])

    # Calculate averages in rings
    total_commuters = all_data.groupby(['smoothed_dist_to_center', 'City'])['B08303e1'].sum()
    total_commute_time = all_data.groupby(['smoothed_dist_to_center', 'City'])['total_commute_time'].sum()
    result = pd.DataFrame(total_commute_time.divide(total_commuters))
    result.reset_index(inplace = True)
    result = result.rename(columns = {0:'avg_commute_time'})
    result = result.loc[result['smoothed_dist_to_center'] <= maximum]

    plot = (ggplot(result, aes(x = 'smoothed_dist_to_center', y = 'avg_commute_time', fill = 'City'))
                   + geom_col(position = 'dodge')
                   + facet_wrap('~City')
                   + theme_bw()
                   + labs(title = 'Commute Times by Distance from the City Center in Austin, Dallas, and Houston',
                          x = 'Distance from City Center (Miles)', y = 'Average Commute Time (Minutes)'))
    plot.save('Figures/Suburbs/travel_times.svg', width = 12, height = 10)

    total_workers = all_data.groupby(['smoothed_dist_to_center', 'City'])['total_workers'].sum()
    total_local_workers = all_data.groupby(['smoothed_dist_to_center', 'City'])['local_workers'].sum()
    workers_pct = 100*pd.DataFrame(total_local_workers.divide(total_workers))
    workers_pct.reset_index(inplace = True)
    workers_pct = workers_pct.rename(columns = {0:'local_workers_pct'})
    workers_pct = workers_pct.loc[workers_pct['smoothed_dist_to_center'] <= maximum]

    plot = (ggplot(workers_pct, aes(x = 'smoothed_dist_to_center', y = 'local_workers_pct', fill = 'City'))
                   + geom_col(position = 'dodge')
                   + facet_wrap('~City')
                   + theme_bw()
                   + labs(title = 'Percent of Workers Working in Place of Residence by Distance from City Center',
                          x = 'Distance from City Center (Miles)', y = 'Percent of Workers Working in Place of Residence'))
    plot.save('Figures/Suburbs/local_workers.svg', width = 12, height = 10)
