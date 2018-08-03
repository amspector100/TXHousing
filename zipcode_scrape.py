import pandas as pd
import geopandas as gpd
import googlemaps
import time
import shapely
from tqdm import tqdm

def correct_dallas_permit_data(address_cache_path = None, dpm_save_path = None):
    """

    Corrects and caches dallas construction permit data. This function worked the first time but is a untested
    since them to be honest. I would not recommend running this again, as it requires getting a google maps account
    and that's been getting harder over time. Instead pull the geocoded addresses from
    https://github.com/amspector100/TXHousing/tree/reorganization/shared_data

    :param address_cache_path: The path at which to cache the geocoded addresses.
    :param dpm_save_path: The path at which to cache the entirely correct permit data.
    :return: None

    """

    import googlemaps
    import shapely
    from tqdm import tqdm

    # Retrieve api key - to get one, follow the directions listed at https://github.com/googlemaps/google-maps-services-python
    api_key_path = "C:/Users/amspe/Desktop/gmaps_apikey.txt"
    with open(api_key_path, 'r') as file:
        api_key = file.read()

    # Instantiate gmaps client
    gmaps = googlemaps.Client(key=api_key)

    def get_latlong(address):
        """

        Get lat and long of an address.

        :param address: String, address. Separated by commas. Example: '1600 Amphitheatre Parkway, Mountain View, CA'
        :return: list(lat, long).
        """
        geocode_result = gmaps.geocode(address)
        return [geocode_result[0]['geometry']['location']['lat'], geocode_result[0]['geometry']['location']['lng']]

    def geocode_dallas_permit_data(permit_types = ['Building (BU) Single Family  New Construction', 'Building (BU) Multi Family  New Construction']):
        """

        When latitudes and longitudes are repeated exactly in the permit data, we assume flag them as possibly incorrect.
        Then, we use the google maps api to re-geocode them. (We cannot geocode every permit as google maps has use limits).
        Note that after the geocode data has been scraped, this will save the geocode data under the save path repeatedly, so that even if
        an error is thrown down the line, you do not need to re-geocode the bad data points.

        :param permit_types: Subset data to only include permit types. Note a lot of multifamily units have repeat geocodes
        (which makes sense) so you run the risk of using too many geocode requests if you put in multifamily new construction repeatedly.
        :param save_path: Save to this path. Should be a shapefile path.
        :return: Corrected dallas permit data.
        """

        bigtime0 = time.time()

        # Get permit data.
        permit_data = process_dallas_permit_data(permit_types = permit_types)

        # Subset based on duplicate lat/long coordinates. This is the part we will correct. Save index.
        subset = permit_data.loc[permit_data['geometry'].apply(lambda x: x.coords[:][0]).duplicated()]

        # Only consider parts with valid addresses. Can't do anything about parts with invalid addresses.
        subset = subset.loc[subset['Street Address'].notnull()]

        # Get unique street addresses. This will reduce the number of necessary geocoding calls. Should be around length 900 for
        # new sf and mf homes.
        addresses = subset['Street Address'].unique().tolist()
        address_dic = {}

        print('Starting to geocode at time {}'.format(time.time() - bigtime0))

        # Loop through. This is a slow method but it is INTENDED to be slow because you can only request 50 geocodes per second.
        for address in tqdm(addresses):

            time0 = time.time()

            # Format address by converting to string (just in case) and adding Dallas, TX info
            address = str(address)
            if 'Dallas, TX' not in address:
                long_address = address + ', Dallas, TX'

            try:
                address_dic[address] = get_latlong(long_address)
            except:
                print('{}, with long_address {}, is a bad address - continuing but you might have to another loop.'.format(address, long_address))

            # Wait if process is working too fast to avoid ratetimeerrors.
            if time.time() - time0 > 0.02:
                time.sleep(0.02 - time.time() + time0)

        data = pd.Series(address_dic)
        data.to_csv(address_cache_path)


    def finish_geocoding():
        """

        Finishes geocoding based on the address dictionary copied into the cached_address_path.txt script.
        Then saves the modified geodata to the dpm_save_path.

        :return: None

        """

        address_dictionary = pd.read_csv(address_cache_path).to_dict()

        # Reget subset
        permit_data = process_dallas_permit_data(permit_types=['Building (BU) Single Family  New Construction', 'Building (BU) Multi Family  New Construction'])
        subset = permit_data.loc[permit_data['geometry'].apply(lambda x: x.coords[:][0]).duplicated()]
        subset = subset.loc[(subset['Street Address'].notnull()) & (subset['Street Address'].isin([key for key in address_dictionary]))]

        # Map to new addresses (tested, this works
        subset['coords'] = subset['Street Address'].map(address_dictionary)
        subset['geometry'] = subset['coords'].apply(lambda x: shapely.geometry.point.Point(x[::-1]))
        subset.drop('coords', axis='columns', inplace=True)


        permit_data.loc[subset.index] = subset

        permit_data.to_file(dpm_save_path)


    # Call this on rows 8825 and 8828 AFTER having read the data in. Google made an error geocoding these rows due to the way
    # it parses addresses (it thinks they're in connecticut).
    def ignore_outliers(rows = [8825, 8828]):
        permit_data = gpd.read_file(dpm_save_path)
        print(permit_data.loc[rows])
        input("Do you want to drop these rows? They should look like they're in Connecticut. If so, type anything - otherwise, stop the program.")
        permit_data.drop(rows, axis = 'index', inplace = True)
        permit_data.to_file(dpm_save_path)

# Utility function: get processed data and rename columns to what their original names were, because the shapefile reader and writer
# in geopandas cuts names off.
def get_corrected_dallas_permit_data(path = dpm_save_path):

    # Old and new column names (new ones are cut off)
    shapefile_columns = ['Permit Num', 'Permit Typ', 'Issued Dat', 'Mapsco', 'Contractor',
       'Value', 'Area', 'Work Descr', 'Land Use', 'Street Add', 'GeoLocatio',
       'Zip Code', 'Year', 'geometry']

    original_columns = ['Permit Number', 'Permit Type', 'Issued Date', 'Mapsco', 'Contractor',
       'Value', 'Area', 'Work Description', 'Land Use', 'Street Address',
       'GeoLocation', 'Zip Code', 'Year', 'geometry']

    dictionary = {}
    for new, old in zip(shapefile_columns, original_columns):
        dictionary[new] = old

    permit_data = gpd.read_file(path)
    permit_data.columns = [dictionary[q] for q in permit_data.columns]
    permit_data.crs = {'init':'epsg:4326'}

    return permit_data



if __name__ == '__main__':
    finish_geocoding()
    ignore_outliers()

