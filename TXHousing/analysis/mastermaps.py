import time
import geopandas as gpd
import folium
from folium import FeatureGroup, LayerControl

# Import other submodules
from .. import utilities
from .. import data_processing
from . import choropleth

def add_zipdata(ziplist, basemap,
                features = ['Views Per Property  (vs CBSA)', 'Hotness Rank ', 'Median Listing Price_sf', 'Median Listing Price_mf'],
                names =  ['Views Per Property (Mean-Centered)', 'Realtor Hotness Rank', 'Median Listing Price for Single Family Homes',
                                                                                 'Median Listing Price for Multifamily Homes']):

    """Adds zipdata to a basemap.


    :param ziplist: A list of zipcodes to add to the basemap
    :param basemap: The basemap to add to.
    :param features: List of features to graph from the realtor datasets.
    :param names: The names of the layers the features will be added on.
    :returns: A ZipBoundaries object with Realtor data attached"""

    zipdata = data_processing.boundaries.ZipBoundaries(ziplist = ziplist)
    zipdata.add_property_data(data_processing.property.realtor_hotness_data, rsuffix = '_hotness')
    zipdata.add_property_data(data_processing.property.realtor_core_inventory_sf, rsuffix = '_sf')
    zipdata.add_property_data(data_processing.property.realtor_core_inventory_mf, rsuffix = '_mf')
    for feature, name in zip(features, names):
        choropleth.continuous_choropleth(zipdata.data, factor = feature, layer_name = name, show = False, basemap = basemap)

    return zipdata


def create_austin_mastermap(save_path = 'Figures/Mastermaps/Austin_Mastermap_v2.html'):

    time0 = time.time()

    basemap = folium.Map([data_processing.zoning.austin_inputs.lat, data_processing.zoning.austin_inputs.long],
                         zoom_start=data_processing.zoning.austin_inputs.zoom)

    # Base zones for municipality. At some point will need to figure out how to plot these.
    regulation_features = ['max_height', 'min_lot', 'far']
    regulation_feature_names = ['Maximum Height', 'Minimum Lot Size (Sq Ft)', 'Maximum Floor to Area Ratio']

    austin_zoning_data = data_processing.zoning.austin_inputs.process_zoning_shapefile(regulation_features = regulation_features)

    # Areas with waived parking requirements
    low_parking_areas = austin_zoning_data.loc[austin_zoning_data['base_zone'].isin(['CBD', 'DMU', 'PD'])]
    choropleth.categorical_choropleth(low_parking_areas, factor = 'base_zone', name = 'Areas with Waived Parking Requirements',
                                      basemap = basemap)

    # Permit and construction data - start by getting permit data
    austin_permit_data = data_processing.permit.process_austin_permit_data(searchfor=['101 single family houses',
                                                                                      '103 two family bldgs',
                                                                                      '104 three & four family bldgs',
                                                                                      '105 five or more family bldgs'],
                                                                           permittypedesc = 'Building Permit', workclass = 'New',
                                                                           earliest=2013)
    def map_permitclass(text):
        return 'SF' if '101' in text else 'MF'

    austin_permit_data['PermitClass'] = austin_permit_data['PermitClass'].apply(map_permitclass)

    # Single family construction - - - - - - - - -
    sfconstruction = austin_permit_data.loc[austin_permit_data['PermitClass'] == 'SF']
    choropleth.make_marker_cluster(sfconstruction, make_centroids = False, fast = True,
                                   name = "Single Family Construction Permits (Marker)", basemap = basemap)
    sfconstruction_grid = utilities.simple.make_point_grid(sfconstruction)
    choropleth.continuous_choropleth(sfconstruction_grid, factor = 'value',
                                      layer_name='Single Family Residential Construction (Choropleth)',
                                      scale_name = 'Number of New Single Family Home Construction Permits in Area (2013-2018)',
                                      show = False, basemap = basemap)
    choropleth.heatmap(sfconstruction, name='Single Family Construction, 2013-2018 (Heatmap)', show=False, radius=13,
                       min_opacity=0.5, max_val=1, basemap = basemap)

    # Multifamily construction - - - - - - - - - -
    mfconstruction = austin_permit_data.loc[austin_permit_data['PermitClass'] == 'MF']
    choropleth.make_marker_cluster(mfconstruction, make_centroids = False, fast = True,
                                   name = "Multifamily Construction Permits (Marker)", basemap = basemap)
    mfconstruction_grid = utilities.simple.make_point_grid(mfconstruction)
    choropleth.continuous_choropleth(mfconstruction_grid, factor = 'value',
                                      layer_name='Multifamily Residential Construction (Choropleth)',
                                      scale_name = 'Number of New Multifamily Home Construction Permits in Area (2013-2018)',
                                      show = False, basemap = basemap)
    choropleth.heatmap(mfconstruction, name='Multifamily Construction, 2013-2018 (Heatmap)', show=False, radius=13, min_opacity=0.5, max_val=1, basemap = basemap)

    # ---------------------------------------------------See print statement---------------------------------------------
    print('Finished permitting layers, took {}. Now creating historic zones markers.'.format(time.time() - time0))

    # Use texas historical sites data for national zones - - - - -
    tx_hd_data = gpd.read_file(data_processing.zoning.tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Austin']
    choropleth.polygon_layer(tx_hd_data, color = 'Green', name = 'National Historic Registry Districts', basemap = basemap)

    # local districts
    signature = '-HD'
    local_hd_data = austin_zoning_data.loc[austin_zoning_data[data_processing.zoning.austin_inputs.feature].str.contains(signature)]
    choropleth.polygon_layer(local_hd_data, color = 'Blue', name = 'Local Historic Districts', basemap = basemap)

    # local landmarks - - - - -
    austin_landmark_data = gpd.read_file(data_processing.zoning.austin_landmark_path)
    choropleth.make_marker_cluster(austin_landmark_data, make_centroids = False, fast = True,
                                   name = 'Local Historic Landmarks', basemap = basemap)

    # ---------------------------------------------------See print statement-------------------------------------------
    print('Finished creating historic layers, took {}. Now creating property layers.'.format(time.time() - time0))

    zipdata = add_zipdata(ziplist = data_processing.boundaries.austin_zips, basemap = basemap)

    # Add dark layer for visualization and layer control
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save(save_path)


def create_dallas_mastermap(save_path = 'Figures/Mastermaps/Dallas_Mastermap_v2.html'):
    """
    Final graph for dallas.

    :param save_path: The path at which to save the final HTML.
    :return: None
    """

    time0 = time.time()

    basemap = folium.Map([data_processing.zoning.dallas_inputs.lat, data_processing.zoning.dallas_inputs.long],
                         zoom_start=data_processing.zoning.dallas_inputs.zoom)

    # Reading initial data
    print('Finished retrieving basemap, took {}. Reading initial data'.format(time.time() - time0))

    base_zones = FeatureGroup('Base Zoning', show = False)

    # Base zones just for Dallas municipality.
    regulation_features = ['height', 'inv_density', 'far']
    regulation_feature_names = ['Maximum Height', 'Minimum Lot Size (Sq Ft)', 'Maximum Floor to Area Ratio']

    dallas_zoning_data = data_processing.zoning.dallas_inputs.process_zoning_shapefile(regulation_features = regulation_features)
    choropleth.categorical_choropleth(dallas_zoning_data, factor = 'broad_zone', colors = ['#D62728', '#17BeCf', '#1f77B4', '#E377C2'], name = 'Base Zoning', basemap = basemap)
    for feature, name in zip(regulation_features, regulation_feature_names):
        print('Working with {}'.format(feature))
        choropleth.continuous_choropleth(dallas_zoning_data.loc[dallas_zoning_data[feature].notnull()],
                                         factor = feature, layer_name=name,
                                         show = False, basemap = basemap)

    choropleth.polygon_layer(dallas_zoning_data.loc[dallas_zoning_data['base_zone'].str.contains('CA-1')],
                             color = 'Blue',
                             name = 'Areas with Weaker Parking Requirements (CA-1 District)',
                             basemap = basemap)

    # Step 2: Historic subdistricts --
    print('In Dallas final graph call, starting to work on historic districts at time {}'.format(time.time() - time0))

    # Start with conservation districts
    cd_districts = dallas_zoning_data.loc[dallas_zoning_data['LONG_ZONE_'].apply(lambda x: x[0:2]) == 'CD']
    choropleth.polygon_layer(cd_districts, color = 'Blue', name = 'Conservation Districts', basemap = basemap)

    # Next do historic overlays
    hist_overlays = gpd.read_file(data_processing.zoning.dallas_historic_overlay_path).to_crs({'init': 'epsg:4326'})
    choropleth.polygon_layer(hist_overlays, color = 'Purple', name = 'Historic Overlays', basemap = basemap)

    # Now do historic subdistricts
    historic_subdistricts = gpd.read_file(data_processing.zoning.dallas_historic_subdistricts_path).to_crs({'init': 'epsg:4326'})
    choropleth.polygon_layer(historic_subdistricts, color = 'Red', name = 'Historic Subdistricts', basemap = basemap)

    # Now do national historic zones
    tx_hd_data = gpd.read_file(data_processing.zoning.tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Dallas']
    choropleth.polygon_layer(tx_hd_data, name = 'National Historic Districts', color = 'Green', basemap = basemap)

    # Step 3: Property/regulatory layers -----------------------------------------------------------------------------
    print('In Dallas final graph call, starting to work on property layers at time {}'.format(time.time() - time0))

    zipdata = add_zipdata(ziplist = data_processing.boundaries.dallas_zips, basemap = basemap)

    # Construction
    all_construction = data_processing.permit.get_corrected_dallas_permit_data()

    # Construction -- single family - - - - - -

    # Marker
    sfconstruction = all_construction.loc[all_construction['Permit Type'] == 'Building (BU) Single Family  New Construction']
    choropleth.make_marker_cluster(sfconstruction, make_centroids = False, fast = True,
                                   name = "Single Family Construction Permits (Marker)", basemap = basemap)
    sfconstruction_grid = utilities.simple.make_point_grid(sfconstruction)
    choropleth.continuous_choropleth(sfconstruction_grid, factor = 'value',
                                      layer_name='Single Family Residential Construction (Choropleth)',
                                      scale_name = 'Number of New Single Family Home Construction Permits in Area (2011-2016)',
                                      show = False, basemap = basemap)
    choropleth.heatmap(sfconstruction, name='Single Family Construction, 2011-2016 (Heatmap)', show=False, radius=13, min_opacity=0.5, max_val=1, basemap = basemap)


    # Construction -- Multifamily - - - - - -
    mfconstruction = all_construction.loc[all_construction['Permit Type'] == 'Building (BU) Multi Family  New Construction']
    choropleth.make_marker_cluster(mfconstruction, make_centroids=False, fast=True,
                                   name = "Multifamily Construction Permits (Marker)", basemap = basemap)
    mfconstruction_grid = utilities.simple.make_point_grid(mfconstruction)
    choropleth.continuous_choropleth(mfconstruction_grid, factor='value',
                                    layer_name='Multifamily Residential Construction (Choropleth)',
                                    scale_name='Number of new Multifamily Construction Permits in Area (2011-2016)',
                                    show=False, basemap = basemap)
    # Heatmap
    choropleth.heatmap(mfconstruction, name='Multifamily Construction, 2011-2016 (Heatmap)', radius=13, min_opacity=0.5, basemap = basemap)


    # Add dark basemap
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save(save_path)


def create_houston_mastermap(save_path = 'Figures/Mastermaps/Houston_Mastermap_v2.html'):
    """Final Houston Mastermap

    :param save_path: The path to save the html to
    :return: None

    """

    basemap = folium.Map([data_processing.zoning.houston_inputs.lat, data_processing.zoning.houston_inputs.long],
                         zoom_start=data_processing.zoning.houston_inputs.zoom)

    # Historic districts -------------------------------------------------------------------

    # National
    tx_hd_data = gpd.read_file(data_processing.zoning.tx_hd_path)
    tx_hd_data = tx_hd_data.loc[tx_hd_data['CITY'] == 'Houston']
    choropleth.polygon_layer(tx_hd_data, name = 'National Historic Districts', color = 'Green', basemap = basemap)

    # Local historic districts
    local_hd_data = gpd.read_file(data_processing.zoning.houston_historic_districts_path).to_crs({'init':'epsg:4326'})
    choropleth.polygon_layer(local_hd_data, name = 'Local Historic Districts', color = 'Blue', basemap = basemap)

    # Local historic landmarks
    local_landmarks_data = gpd.read_file(data_processing.zoning.houston_historic_landmarks_path).to_crs({'init':'epsg:4326'})
    choropleth.make_marker_cluster(local_landmarks_data, make_centroids = False, fast = True,
                                   name = 'Local Historic Landmarks', basemap = basemap)

    # Construction permit data ----------------------------------------------------------------------
    print('Processing Houston permit data')
    houston_permit_data = data_processing.permit.process_houston_permit_data(searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE',
                                                                                           'NEW AP', 'NEW HI-'],
                                                                              searchin = ['PROJ_DESC'],
                                                                              earliest = 2013, latest = None)

    # Subset to only include approved permits and nonempty geometries
    houston_permit_data = houston_permit_data.loc[houston_permit_data['Approval'] == 1.0]
    houston_permit_data = utilities.simple.process_points(houston_permit_data)

    # SF construction - marker cluster, choropleth, and heatmap
    sfconstruction = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW S.F.', 'NEW SF', 'NEW TOWNHOUSE', 'NEW SINGLE']))]
    choropleth.make_marker_cluster(sfconstruction, make_centroids = False, fast = True,
                                   name = "Single Family Construction Permits (Marker)", basemap = basemap)
    sfconstruction_grid = utilities.simple.make_point_grid(sfconstruction)
    choropleth.continuous_choropleth(sfconstruction_grid, factor = 'value',
                                      layer_name='Single Family Residential Construction (Choropleth)',
                                      scale_name = 'Number of New Single Family Home Construction Permits in Area (2013-2018)',
                                      show = True, basemap = basemap)

    choropleth.heatmap(sfconstruction, name='Single Family Construction, 2013-2018 (Heatmap)', show=False, radius=13,
                       min_opacity=0.5, max_val=1, basemap = basemap)

    # MF construction
    mfconstruction = houston_permit_data.loc[houston_permit_data['PROJ_DESC'].str.contains('|'.join(['NEW AP', 'NEW HI-']))]
    choropleth.make_marker_cluster(mfconstruction, make_centroids=False, fast=True,
                                   name = "Multifamily Construction Permits (Marker)", basemap = basemap)
    mfconstruction_grid = utilities.simple.make_point_grid(mfconstruction)
    choropleth.continuous_choropleth(mfconstruction_grid, factor='value',
                                    layer_name='Multifamily Residential Construction (Choropleth)',
                                    scale_name='Number of new Multifamily Construction Permits in Area (2013-2018)',
                                    show=False, basemap = basemap)
    # Heatmap
    choropleth.heatmap(mfconstruction, name='Multifamily Construction, 2013-2018 (Heatmap)', radius=13,
                       min_opacity=0.5, basemap = basemap)

    # Add dark layer for visualization, layer control, then save
    folium.TileLayer('cartodbdark_matter').add_to(basemap)
    LayerControl().add_to(basemap)
    basemap.save(save_path)
