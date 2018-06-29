# Some global paths ---------------------------------------------------------------------------------------------------
nbhd_boundaries_path = "data\Zillow Data\ZillowNeighborhoods-TX\ZillowNeighborhoods-TX.shp"
zip_boundaries_path = "data\cb_2017_us_zcta510_500k\cb_2017_us_zcta510_500k.shp"

# Class for base zoning inputs ----------------------------------------------------------------------------------------
class zoning_inputs():

    def __init__(self, path, feature, separator, proj4string, base_zones, lat = 0, long = 0, zoom = 10,
                 title = '', xlims = [None, None], ylims = [None, None], crs = None):
        self.path = path
        self.feature = feature
        self.separator = separator
        self.proj4string = proj4string
        self.base_zones = base_zones
        self.lat = lat
        self.long = long
        self.zoom = zoom
        self.title = title
        self.xlims = xlims
        self.ylims = ylims
        self.crs = crs

# Downloaded from http://data-nctcoggis.opendata.arcgis.com/datasets/2015-land-use
north_texas_inputs = zoning_inputs(path = 'data/Zoning Shapefiles/2015_North_Texas_Land_Use/2015_Land_Use.shp',
                                   feature = 'CATEGORY',
                                   separator = 'no_separator',
                                   proj4string = 'EPSG:4326',
                                   base_zones = { 'Single Family':['Single family'],
                                               'Multifamily':['Mixed use', 'Multi-family'],
                                               'Other Residential':['Residential acreage'],
                                               'Other':[]},
                                   lat = 32.7767,
                                   long = -96.7970,
                                   zoom = 10,
                                   title = 'Base Zones in and around Dallas, Texas')
                                   #xlims = [-97.3, -96.1],
                                   #ylims = [32.25, 33.25])
# Downloaded from https://data-nctcoggis.opendata.arcgis.com/datasets/81916863ab394786ab5caaa731f5ac36_4?geometry=-104.084%2C30.991%2C-83.342%2C34.229
nt_highways = zoning_inputs(path = "data/Zoning Shapefiles/NT_Highways_2017/Highways_2017.shp",
                            feature = 'unsure',
                            separator = 'no_separator',
                            proj4string = 'EPSG:4326',
                            base_zones = {},
                            long = north_texas_inputs.long,
                            lat = north_texas_inputs.lat,
                            zoom = north_texas_inputs.zoom)

# Downloaded from https://data-nctcoggis.opendata.arcgis.com/datasets/counties-
nt_counties = zoning_inputs(path = "data/Zoning Shapefiles/nt_counties/Counties_.shp",
                            feature = '',
                            separator = '',
                            proj4string = 'EPSG:4326',
                            base_zones = {},
                            long = north_texas_inputs.long,
                            lat = north_texas_inputs.lat,
                            zoom = north_texas_inputs.zoom)



# Downloaded from https://data.austintexas.gov/Locations-and-Maps/Zoning/5rzy-nm5e
austin_inputs = zoning_inputs(path = "data/Zoning Shapefiles/austin_zoning/geo_export_571668ee-52f1-4ac9-a4e0-3b8bb348eae7.shp",
                              feature = 'zoning_zty',
                              separator = '-',
                              proj4string = 'EPSG:4326',
                              base_zones = {'Single Family':['SF'], 'Multifamily':['MF'],
                                            'Other Residential':['MH', 'RR', 'LA'], 'Other':[]},
                              lat = 30.267,
                              long = -97.743,
                              zoom = 9,
                              title = 'Base Zones in Austin, Texas',
                              crs = {'init': 'epsg:4326'})

# From dallas's open data site
dallas_inputs = zoning_inputs(path = "data/Zoning Shapefiles/DallasBaseZoning/BaseZoning.shp",
                              feature = 'ZONE_DIST',
                              separator = '-',
                              proj4string = 'EPSG:2276',
                              base_zones = {'Single Family':['A(A)', 'D(A)', 'TH', 'R'],
                                            'Multifamily': ['CH', 'MF', 'MU'],
                                            'Other Residential':['MH(A)'],
                                            'Other':[]},
                              lat = north_texas_inputs.lat,
                              long = north_texas_inputs.long,
                              zoom = 9,
                              title = 'Base Zones in Dallas, Texas')



# City class  --------------------------------------------------------------------------------------------------------
class City:
    def __init__(self, name, long, lat, zoom, xlims = [None, None], ylims = [None, None]):
        self.name = name
        self.long = long
        self.lat = lat
        self.zoom = zoom
        self.xlims = xlims
        self.ylims = ylims

# Cities (needed for defining Demand_Input)
austin = City('Austin', lat =30.267, long =-97.743, zoom=11, xlims = [29.25, 31.75], ylims = [96.5, 99])
dallas = City('Dallas', lat = 32.7767, long =-96.7970, zoom = 11, xlims = [-97.3, -96.1], ylims = [32.25, 33.25])
houston = City('Houston', lat = 29.7604, long  = -95.3698, zoom = 11, xlims = [28, 31], ylims = [-93, -97])

# Demand input class ------------------------------------------------------------------------------------------------
class Demand_Input:

    def __init__(self, path, source,
                 feature = 'Value',
                 geo_filter_classes=[austin, dallas, houston],
                 geography = 'zip',
                 **kwargs):
        """
        :param path: Path of the dataset, assumed csv
        :param source: Source of the data. Can either be "Zillow" or "Realtor"
        :param feature: Name of the feature. In the case of Realtor data, this must be the column of the data.
        :param geo_filter_classes: A list of City classes which are used to subset the data.
        :param geography: The level of detail of the dataset, i.e. "zip" or "neighborhood."
        :param geo_filter: A kwarg. This is the geography level (i.e. city, state, county) by which to filter data,
        and it should be a column of the dataset. Default depends on the source.
        :param index_col: A kwarg. The column to use as the index. Default depends on the source.
        :param name: an alternate name for graphical display of the feature. Defaults to the feature string.
        """


        self.path = path
        self.source = source
        self.feature = feature
        self.geo_filter_classes = geo_filter_classes
        self.geography = geography

        # Depending on the source, use different defaults
        try:
            self.geo_filter == kwargs['geo_filter']
        except:
            if source == 'Zillow':
                self.geo_filter = 'City'
            elif source == 'Realtor':
                self.geo_filter = 'ZipName'

        try:
            self.index_col = kwargs['index_col']
        except:
            if source == 'Zillow':
                self.index_col = 'RegionName'
            elif source == 'Realtor':
                self.index_col = 'ZipCode'


        # Get filter values and name from kwargs
        try:
            self.geo_filter_values = kwargs['geo_filter_values']
        except:
            self.geo_filter_values = [city.name for city in self.geo_filter_classes]

        try:
            self.name = kwargs['name']
        except:
            self.name = feature

# Zillow inputs
sfhomes_nbhd = Demand_Input(path = 'data/Zillow Data/sfhomes_neighborhood.csv',
                        source = 'Zillow',
                        feature = 'Single Family Homes Median Value',
                        geography = 'Neighborhood')

inventoryraw_zip = Demand_Input(path = "data/Zillow Data/inventoryraw_zip.csv",
                            source = 'Zillow',
                            feature = 'Raw Inventory',
                            geography = 'Zip')

medpricecuts_zip = Demand_Input(path = "data/Zillow Data/medpricecut_zip.csv",
                            source = 'Zillow')

allhomeprices_zip = Demand_Input(path = 'data/Zillow Data/allhomeprice_zip.csv',
                             source = 'Zillow')

# One note is that to avoid a rather annoying index error I had to change the name of "Downtown" region in Austin to
# "downtown" (without the capital D).
allhomeprices_nbhd = Demand_Input(path = "data/Zillow Data/allhomeprice_neighborhood.csv",
                                  source = 'Zillow',
                                  geography = 'Neighborhood')

# Realtor inputs
realtor_med_dom_cbsa = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Median DOM (vs CBSA)',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_med_dom_yoy = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Median DOM Y/Y',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_med_vpp_cbsa = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Views Per Property  (vs CBSA)',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])


realtor_med_vpp_yoy = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Views Per Property Y/Y',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_hotness_cbsa = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Hotness Rank Within CBSA',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_hotness_yoy = Demand_Input(path = "data/RDC_MarketHotness_Monthly_Zip.csv",
                               source = 'Realtor',
                               feature='Hotness Rank Y/Y',
                               geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

# Core realtor data
realtor_avg_price = Demand_Input(path = 'data/RDC_InventoryCoreMetrics_Zip.csv',
                                 source = 'Realtor',
                                 feature = 'Avg Listing Price',
                                 geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'])

realtor_avg_sf_price = Demand_Input(path = "data/RDC_InventoryCoreMetrics_Zip_sfh.csv",
                                 source = 'Realtor',
                                 feature = 'Avg Listing Price',
                                 geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'],
                                 name='Avg Listing Price, SF Homes')


realtor_avg_cth_price = Demand_Input(path = "data/RDC_InventoryCoreMetrics_Zip_cth.csv",
                                 source = 'Realtor',
                                 feature = 'Avg Listing Price',
                                 geo_filter_values=['Dallas, TX', 'Austin, TX', 'Houston, TX'],
                                 name='Avg Listing Price, Condo / Townhomes')

# Some more specific inputs --------------------------------------------------------------------------------------------

# Austin -------------------------------------------------------------------------------------------------------------

# Parcel data
austin_parcel_path = "data/Zoning Shapefiles/Austin Land Database 2016/geo_export_813e97e4-7fde-4e3a-81b3-7ca9e8a89bd0.shp"

# Zip codes in austin
austin_zips = [73301, 73344, 78704, 78705, 78708, 78713, 78714, 78715, 78701, 78702, 78703, 78709, 78710, 78711, 78712,
78716, 78717, 78718, 78719, 78720, 78721, 78722, 78723, 78728, 78729, 78730, 78731, 78734, 78724, 78725, 78726, 78727,
78732, 78733, 78735, 78736, 78739, 78741, 78745, 78746, 78749, 78750, 78751, 78752, 78755, 78756, 78757, 78758, 78759,
78737, 78738, 78742, 78744, 78747, 78748, 78753, 78754, 78760, 78761, 78762, 78763, 78769, 78772, 78773, 78780, 78781,
78799, 78764, 78765, 78766, 78767, 78768, 78774, 78778, 78779, 78783, 78785, 78789]
austin_zips = [str(zip) for zip in austin_zips]

# Regulation features and fill values. The fill values for: max_units_lot, max_units_acre, max_height, and far, are in question.
austin_regulation_types = {"min_lot":0, "min_lot_wid":0, "max_units_lot":10, "max_units_acre":10, "max_build_cov":100, "max_imperv":100,
                           "max_height":200, "front_yard":0, "street_yard":0, "interior_yard":0, "rear_yard":0, "far":15,
                           "min_site_area_unit":0, "min_site_area_build":0}


austin_zip_features_dic = {'Median DOM (vs CBSA)':['Median Days on Market', 'Median Days on Market, Data from Realtor'],
                    'far':['Floor to Area Ratio', 'Floor to Area Ratio, Data from City of Austin'],
                    'min_lot':['Minimum Lot Size', 'Minimum Lot Size (Square Feet)'],
                    'max_build_cov':['Maximum Building Coverage', 'Maximum Building Coverage (%)'],
                   'sf_avg_listing': ['Single Family Homes Price', 'Average Single Family Home Listing Price (Realtor)'],
                   'mf_avg_listing': ['Multifamily Homes Price', 'Average Multifamily Home Listing Price (Realtor)']}

austin_processed_parcel_data_path = 'data/Processed Austin Parcel Data/processed_parcel_data.shp'
austin_downzone_path = 'data/Zoning Shapefiles/Austin Downzoning/austin_downzoning.shp'
austin_construction_path = 'data/Zoning Shapefiles/Austin Construction/austin_construction.shp'
austin_permit_path = 'data/austin_construction_permits.csv'

austin_regulations_path = "data/austin zoning standards.csv"

# Dallas -------------------------------------------------------------------------------------------------------------

dallas_zips = [75203, 75204, 75205, 75208, 75209, 75210, 75211, 75212, 75214, 75201, 75202, 75206, 75207, 75215, 75216,
75217, 75218, 75222, 75223, 75224, 75225, 75230, 75231, 75232, 75233, 75236, 75237, 75238, 75240, 75246, 75251, 75252,
75253, 75219, 75220, 75221, 75226, 75227, 75228, 75229, 75234, 75235, 75241, 75242, 75243, 75244, 75247, 75248, 75249,
75250, 75254, 75260, 75261, 75265, 75266, 75267, 75270, 75275, 75284, 75285, 75262, 75263, 75264, 75277, 75283, 75287,
75301, 75315, 75320, 75326, 75355, 75356, 75357, 75358, 75370, 75371, 75372, 75373, 75381, 75382, 75389, 75390, 75393,
75394, 75395, 75303, 75312, 75313, 75336, 75339, 75342, 75354, 75359, 75360, 75367, 75368, 75374, 75376, 75378, 75379,
75380, 75391, 75392, 75397, 75398]

dallas_zips = [str(zip) for zip in dallas_zips]

# The fill values for: far and height are in question
dallas_regulation_types = {'front':0, 'side':0, 'inv_density':0, 'far':15, 'height':200, 'coverage':100}

dallas_zip_features_dic = {'Median DOM (vs CBSA)':['Median Days on Market', 'Median Days on Market, Data from Realtor'],
                           'far': ['Floor to Area Ratio', 'Floor to Area Ratio, Data from City of Dallas'],
                           'inv_density': ['Minimum Dwelling Size', 'Minimum Dwelling Size (Square Feet)'],
                           'coverage': ['Maximum Building Coverage', 'Maximum Building Coverage (%)'],
                           'sf_avg_listing': ['Single Family Homes Price', 'Average Single Family Home Listing Price (Realtor)'],
                           'mf_avg_listing': ['Multifamily Homes Price', 'Average Multifamily Home Listing Price (Realtor)']}

dallas_parcel_data_path_2013 = "data/Zoning Shapefiles/Dallas County Parcels 2013/geo_export_9b090abf-d5d9-4c74-a6be-4486e75ee147.shp"
dallas_parcel_data_path_2016 = "data/Zoning Shapefiles/Dallas County Parcels 2016/geo_export_bd65f212-41ce-4166-804a-5dc5ef85ee84.shp"

dallas_processed_parcel_data_path = 'data/Processed Dallas Parcel Data/processed_parcel_data.shp' # This doesn't work for some reason
dallas_downdev_path = 'data/Zoning Shapefiles/Dallas Downdevelopments/dallas_downdevelopments.shp'
dallas_construction_path = 'data/Zoning Shapefiles/Dallas Construction/dallas_construction.shp'
dallas_permit_path = 'data/dallas_construction_permits.csv'

dpm_save_path = "C:/Users/amspe/Documents/R/MI2018/TXHousing/data/Zoning Shapefiles/Dallas Corrected Permits/dallas_permits_corrected.shp"

dallas_regulations_path = 'data/dallas zoning standards.csv'

renovation_types_dallas = ['Building (BU) Commercial  Renovation', 'Building (BU) Single Family  Renovation',
                           'Building (BU) Multi Family  Renovation', 'Flammable Liquid (FL) Commercial  Renovation',
                           'Swimming Pool (SW) Single Family  Renovation', 'Tree Removal (TR) Commercial  Renovation',
                           'Swimming Pool (SW) Multi Family  Renovation', 'Tent (TE) Commercial  Renovation',
                           'Excavation (EX) Commercial  Renovation',
                           'Grading and Paving (GP) Single Family  Renovation',
                           'Grading and Paving (GP) Commercial  Renovation', 'Electrical (EW) Commercial  Renovation',
                           'Plumbing (PL) Commercial  Renovation', 'Fence (FE) Single Family  Renovation',
                           'Liquor License (LL) Commercial  Renovation', 'Swimming Pool (SW) Commercial  Renovation',
                           'Tree Removal (TR) Multi Family  Renovation', 'Excavation (EX) Single Family  Renovation',
                           'Tent (TE) Single Family  Renovation', 'Security System (SE) Commercial  Renovation',
                           'Plumbing (PL) Single Family  Renovation', 'Tree Removal (TR) Single Family  Renovation',
                           'Electrical (EL) Commercial  Renovation', 'Fence (FE) Multi Family  Renovation',
                           'Fire Alarm (FA) Commercial  Renovation', 'Fence (FE) Commercial  Renovation',
                           'Grading and Paving (GP) Multi Family  Renovation',
                           'Paving (Sidewalk, Drive Approaches) (PV) Single Family  Renovation',
                           'Mechanical (ME) Single Family  Renovation', 'Tree Mitigation (MI) Commercial  Renovation',
                           'Barricade (BA) Commercial  Renovation', 'Electrical (EW) Single Family  Renovation',
                           'Mechanical (ME) Commercial  Renovation', 'Lawn Sprinkler (LS) Single Family  Renovation',
                           'Electrical (EL) Single Family  Renovation',
                           'Fire Sprinkler (Minor Work) (FS) Commercial  Renovation',
                           'Electrical (EL) Multi Family  Renovation', 'Mechanical (ME) Multi Family  Renovation',
                           'Plumbing (PL) Multi Family  Renovation', 'Lawn Sprinkler (LS) Commercial  Renovation',
                           'Excavation (EX) Multi Family  Renovation', 'Flammable Liquid (FL) Multi Family  Renovation',
                           'Fire Sprinkler (Minor Work) (FS) Single Family  Renovation',
                           'Electrical (EW) Multi Family  Renovation',
                           'Fire Sprinkler (Minor Work) (FS) Multi Family  Renovation',
                           'Medical Gas (MG) Commercial  Renovation', 'Backflow (BF) Commercial  Renovation',
                           'Paving (Sidewalk, Drive Approaches) (PV) Commercial  Renovation',
                           'Barricade (BA) Multi Family  Renovation', 'Tent (TE) Multi Family  Renovation']


# Misc -------------------------------------------------------------------------------------------------------------

# See metadata at https://www2.census.gov/geo/tiger/TIGER_DP/2016ACS/Metadata/BG_METADATA_2016.txt
texas_blocks_path = "data/ACS_2016_5YR_BG_48_TEXAS.gdb"
texas_places_path = "data/cb_2017_48_place_500k/cb_2017_48_place_500k.shp"

# And colors - divide by 256 for consistency with matplotlib
import numpy as np
white = np.array((256, 256, 256), dtype = int)/256
blue = np.array((90, 90, 250))/256
green = np.array((88, 170, 168), dtype = int)/256
red = np.array((250, 90, 90))/256
yellow = np.array((240, 240, 90))/256
purple = np.array((18, 15, 67), dtype = int)/256