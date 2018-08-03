import os
import time
import pandas as pd
import geopandas as gpd
import warnings
import shapely.geometry


# Filepaths
austin_permit_path = 'data/austin_construction_permits.csv'

dallas_permit_path = 'data/dallas_construction_permits.csv'
dpm_save_path = "data/Zoning Shapefiles/Dallas Corrected Permits/dallas_permits_corrected.shp"

houston_structural_permits_path = "data/Houston_Structural_Permits/Permits_wm_Structural.shp"
houston_demolition_permits_path = "data/Houston_Demolition_Permits/Demolition_ILMS_Code_SD.shp"
houston_permit_statuses_path = 'shared_data/houston_permit_statuses.csv'
backup_houston_permit_statuses_path = 'shared_data/houston_permit_statuses_backup.csv'

# Austin -----------------------------------

def process_austin_permit_data(searchfor, permittypedesc=None, workclass=None,
                               earliest=None, latest=None):
    """
    :param searchfor: List of strings to search for in the 'PermitClass' and 'Description' columns. Often it's worth
        explicitly passing in permit classes, i.e. searchfor = ['101 single family houses'].
    :param permittypedesc: The permittypedesc to match. Ex: "Building Permit."
    :param workclass: Workclass to match. Ex: "New"
    :param earliest: Earliest date to consider (inclusive). Data runs from 1971-2018.
    :param latest: Latest date to consider (inclusive). Data runs from 1971-2018.
    :return: GeoDataFrame of subsetted permit data.
    """

    time0 = time.time()
    print('Reading Austin permit data')
    permit_data = pd.read_csv(austin_permit_path)
    print('Finished reading Austin permit data, took {}. Now subsetting.'.format(time.time() - time0))

    # Initial subset
    construction = permit_data
    if permittypedesc is not None:
        construction = construction.loc[(permit_data['PermitTypeDesc'] == permittypedesc)]
    if workclass is not None:
        construction = construction.loc[(permit_data['WorkClass'] == workclass)]

    # Get right timeframe
    if earliest is not None:
        construction = construction.loc[construction['CalendarYearIssued'] >= earliest, :]

    if latest is not None:
        construction = construction.loc[construction['CalendarYearIssued'] <= latest, :]

    def to_lowercase(object):
        return str(object).lower()

    # Now subset to just include entries where one of the 'searchfor' strings is in either the description or permit
    #  class column. Can use this to search for either multifamily or single family, or anything specific really.
    # Best to search based on permit class, see below.
    searchfor = [str(x).lower() for x in searchfor]
    construction.loc[:, 'PermitClass'] = construction.loc[:, 'PermitClass'].apply(to_lowercase)
    construction.loc[:, 'Description'] = construction.loc[:, 'Description'].apply(to_lowercase)
    construction = construction.loc[(construction['PermitClass'].str.contains('|'.join(searchfor))) |
                                    (construction['Description'].str.contains('|'.join(searchfor))), :]

    def return_point(series):
        # series[0] is long, series[1] is lat
        return shapely.geometry.point.Point(series[0], series[1])

    # Now process geometry
    points = construction[['Longitude', 'Latitude']].apply(return_point, axis=1)
    construction = gpd.GeoDataFrame(data=construction, geometry=points)
    construction.crs = {'init': 'epsg:4326'}

    print('Finished processing Austin permit data, took {}.'.format(time.time() - time0))
    return construction

# Dallas -----------------------------------------------

def process_dallas_permit_data(permit_types, earliest=None, latest=None):
    """
    Initially process raw dallas permit data. However, the raw permit data has some inaccuracies, and
    it's best to simply read the corrected data from dpm_save_path.

    :param permit_types: List of permit types to filter for. Will only consider rows where the permit type is one of
        these permit types.
    :param earliest: Earliest date to consider. Data runs from 2011-2016. Defaults to None.
    :param latest: Latest date to consider. Data runs from 2011-2016. Defaults to None.
    :raises: UserWarning; some permit data is incorrect. Construction data has been corrected, and the corrected data is
        stored at 'dpm_save_path' as listed in the inputs.py file.
    :return: GeoDataFrame of subsetted permit data
    """

    warnings.warn("""Some permit data is incorrect. Construction permit data has been corrected and the corrected data 
     should be stored at "dpm_save_path" file. Instead of running this function, you should simply read that data using
     the get_corrected_dallas_permit_data function.""")

    permit_data = pd.read_csv(dallas_permit_path)

    # Get date
    def get_year(text):
        date = text.split(' ')[0]
        year = date.split('/')[-1]
        return int(year)

    permit_data['Year'] = permit_data['Issued Date'].apply(get_year)

    # subset by year and permit type, then return
    construction = permit_data.loc[permit_data['Permit Type'].isin(permit_types)]
    if earliest is not None:
        construction = construction.loc[construction['Year'] >= earliest]
    if latest is not None:
        construction = construction.loc[construction['Year'] <= latest]

    # Get points
    def return_point(location):
        coords = str(location).split('\n')[-1]
        lat = float(coords[1:].split(', ')[0])
        long = float(coords.split(', ')[-1][:-1])
        return shapely.geometry.point.Point(long, lat)

    # Geocoding fails for 5000/220,000 locations, which have addresses but not lat/long. Could improve this using google API
    # or using a dictionary of addresses:lat/long (because some addresses are repeated) but a bit pressed for time right now.
    # Fails for about 15% of construction sites.
    def is_geocoded(location):
        try:
            return_point(location)
            return True
        except:
            return False

    construction = construction.loc[construction['GeoLocation'].apply(is_geocoded)]
    points = construction['GeoLocation'].apply(return_point)
    construction = gpd.GeoDataFrame(data=construction, geometry=points)

    return construction


def correct_dallas_permit_data(address_cache_path = None, dpm_save_path = None, re_geocode = True):
    """

    Corrects and caches dallas construction permit data. This function worked the first time but is a untested
    since then. I would not recommend running this again, as it requires getting a google maps account and that's been
    getting harder over time. Instead, it's easier to pull the geocoded addresses from
    https://github.com/amspector100/TXHousing/tree/reorganization/shared_data

    :param address_cache_path: The path at which to cache the geocoded addresses.
    :param dpm_save_path: The path at which to cache the entirely correct permit data.
    :param re_geocode: If true, re-pull all the data from the google maps API. Else, just use the cached address path.
    :type re_geocode: Boolean
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

    def finish_geocoding(permit_types = ['Building (BU) Single Family  New Construction', 'Building (BU) Multi Family  New Construction']):

        # Start over from this point to make everything dependent on the cached data and not the geocoding.
        address_dictionary = pd.read_csv(address_cache_path).to_dict()

        # Reget subset
        permit_data = process_dallas_permit_data(permit_types=['Building (BU) Single Family  New Construction', 'Building (BU) Multi Family  New Construction'])
        subset = permit_data.loc[permit_data['geometry'].apply(lambda x: x.coords[:][0]).duplicated()]
        subset = subset.loc[(subset['Street Address'].notnull()) & (subset['Street Address'].isin([key for key in address_dictionary]))]

        # Map to new addresses. There are also some rows which will probably show up as if they're in Connecticut due
        # to address parsing error on Google's end - it will offer you the chance to drop these points.
        subset['coords'] = subset['Street Address'].map(address_dictionary)
        subset['geometry'] = subset['coords'].apply(lambda x: shapely.geometry.point.Point(x[::-1]))
        subset.drop('coords', axis='columns', inplace=True)
        permit_data.loc[subset.index] = subset


        rows = [8825, 8828]
        print(permit_data.loc[rows])
        flag = input("Do you want to drop the rows which were just printed? They look like they're in Connecticut. If so, type 'yes'.")
        if flag.lower() == 'yes':
            permit_data.drop(rows, axis='index', inplace=True)


        permit_data.to_file(dpm_save_path)

    if re_geocode:
        geocode_dallas_permit_data()
    finish_geocoding()

def get_corrected_dallas_permit_data(path = dpm_save_path):
    """Get processed & corrected Dallas permit data.

    :param path: the path to read the data from.
    :return: GDF of Construciton permit data.
    """

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

    try:
        permit_data.columns = [dictionary[q] for q in permit_data.columns]
    except:
        warnings.warn('Corrected dallas permit data column names have been cut off by shapefile editor')

    permit_data.crs = {'init':'epsg:4326'}

    return permit_data


# Houston ---------------------------------------------------------------------

def process_houston_permit_data(searchfor=['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE'],
                                searchin=['PROJ_DESC'],
                                kind='structural',
                                earliest=None, latest=None):
    """ Process houston permit data. Note that houston permit data does not specify whether housing is new or not, this
    requires parsing the descriptions to figure out whether housing is new sf/new mf housing. Thankfully the descriptions
    are formulaic so this is not too hard to do.

    :param searchfor: A list of keywords to search for, i.e. ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE']).
        To subset to only include new construction, each string should start with 'NEW '.
    :param searchin: The columns to search in. Defualts to ['PROJ_DESC'], project description. Note that the function
        will return rows where ANY of the columns specified by 'searchin' contain ANY of the keywords in 'searchfor'.
    :param kind: The kind of permit data to read in. Defaults to 'structural', can either by 'structural'
        or 'demolition'
    :param earliest: Earliest date to consider.
    :param latest: Latest date to consider.
    :return: GeoDataFrame of subsetted permit data.
    """

    # Load data
    if kind == 'structural':
        path = houston_structural_permits_path
    elif kind == 'demolition':
        path = houston_demolition_permits_path
    else:
        raise ValueError(""""kind" arg must either equal "demolition" or "structural", not "{}""".format(kind))

    permit_data = gpd.read_file(path)
    permit_data = permit_data.rename(
        columns={'SITUS_ZIP_': 'Zipcode'})  # Already has 100% complete zipcode information.
    permit_data = permit_data.loc[~permit_data['geometry'].apply(lambda x: x is None)]

    # Ignore repeated plans
    permit_data = permit_data.loc[~permit_data['PROJ_DESC'].str.contains('REPEAT PLAN')]

    # Search for 'searchfor' args in the specified columns
    flags = pd.Series(False, index=permit_data.index)
    for col in searchin:
        flags = flags | permit_data[col].str.contains('|'.join(searchfor))
    subsetted_data = permit_data.loc[flags]

    # Process dates. Note that data runs from 1988 to 2018.
    subsetted_data.loc[:, 'Year'] = subsetted_data['APPLN_DATE'].apply(lambda x: int(x[0:4]))
    if earliest is not None:
        subsetted_data = subsetted_data.loc[subsetted_data['Year'] >= earliest]
    if latest is not None:
        subsetted_data = subsetted_data.loc[subsetted_data['Year'] <= latest]

    # Process whether permit has been approved yet (using scraped data)
    if os.path.exists(houston_permit_statuses_path):
        approvals = pd.read_csv(houston_permit_statuses_path, index_col=0)
        approvals.index = [str(id) for id in approvals.index]
    else:
        warnings.warn("""'Cannot retrieve accurate approval status for houston permit data because it has not been
         scraped or downloaded. Try pulling it from https://github.com/amspector100/TXHousing/tree/reorganization/shared_data
         or alternatively run the scrape_houston_permit_data function.""")

    # Merge
    subsetted_data['Approval'] = subsetted_data['PROJECT_NO'].map(approvals['Approval'])

    # Return
    return subsetted_data

def scrape_houston_permit_data(target_path = houston_permit_statuses_path, kind = 'structural', backup_path = backup_houston_permit_statuses_path):
    """Scrapes Houston permit approval data and writes it to the target_path as well as the backup_path. It will never
    overwrite any pre-existing csv files - it only appends information to csvs. This requires proper installation of the
    headless chrome webdriver and takes a lot of time (~2 hours), so it's probably best to pull the information from
    https://github.com/amspector100/TXHousing/tree/reorganization/shared_data instead of calling this function.

    :param target_path: The path to write the scraped data to.
    :param kind: The kind of permit data to scrape. Can either be 'structural' or 'demolition'.
    :param backup_path: The path to write the scraped data to as a backup in case something happens to the original.
    :returns: None.
    """

    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
    from multiprocessing import Pool
    import sys

    def get_project_statuses(project_numbers):
        """
        Scrapes status of permit projects from Houston website.
        :param project_numbers: List of project numbers (strings).
        :return: A pandas series. The index is the project numbers, the values are the approval status.
        0 = not approved, 1 = approved, Na = no information listed on the website.
        """
        result_dictionary = {}
        counter = 0

        # Open webdriver and initialize
        options = webdriver.ChromeOptions()
        options.add_argument('--headless')
        options.add_argument('--log-level=3')
        browser = webdriver.Chrome(chrome_options=options)
        url = 'https://www.pdinet.pd.houstontx.gov/cohilms/webs/Plan_LookUp.asp'

        for project_number in project_numbers:

            if counter % 100 == 0:
                print('{} process has finished {} ids'.format(__name__, counter))

            # Some very crude error handling
            try:
                # Go to website
                browser.get(url)

                # Search for specific project
                input_bar = browser.find_element_by_name('ProjectNo')
                input_bar.send_keys(project_number)
                input_bar.send_keys(Keys.RETURN)  # Hit return after sending text

                # Get content and return result
                content = browser.find_elements_by_class_name('content')
                content_list = [item.text for item in content]
            except:
                print("Unexpected error for {}".format(project_number), sys.exc_info()[0])
                return result_dictionary

            if 'These plans have NOT been approved for permitting' in content_list:
                result_dictionary[project_number] = 0.0  # Not approved
            elif len(content_list) <= 3:
                result_dictionary[project_number] = float('NaN')  # No data
            else:
                result_dictionary[project_number] = 1.0  # Approved

            counter += 1

        # Quit browser
        browser.quit()
        result_dictionary = pd.Series(result_dictionary)
        return result_dictionary

    # Use multithreading to do the same thing but faster
    def speedy_get_project_statuses(all_project_numbers, number_processes=3):

        time0 = time.time()
        chunksize = round(len(all_project_numbers) / number_processes) + 1

        def chunks(l, n):
            """Yield successive n-sized chunks from l.
            :param l: The list in question
            :param n: the chunksize."""
            for i in range(0, len(l), n):
                yield l[i:i + n]

        all_project_numbers_chunked = list(chunks(all_project_numbers, chunksize))

        p = Pool(number_processes)
        results = p.map(get_project_statuses, all_project_numbers_chunked)
        p.close()
        p.join()
        results = pd.concat(results)
        results.to_csv(target_path, mode='a', header=False)
        print('Took {} seconds for {} project numbers'.format(time.time() - time0, len(all_project_numbers)))
        return results

    if __name__ == '__main__':

        from helpers import process_houston_permit_data

        # Step 1: Get all the unique project numbers
        if kind == 'structural':
            searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE', 'NEW AP', 'NEW HI-'] # All residential buildings
        elif kind == 'demolition':
            searchfor = [''] # Almost everything is residential demolition, and anyway, commercial demos are included in Dallas

        searchin = ['PROJ_DESC']
        kind = 'demolition'
        earliest = 2010
        latest = None

        print('Processing permit data')
        houston_permit_data = process_houston_permit_data(searchfor=searchfor,
                                                          searchin=searchin,
                                                          kind=kind,
                                                          earliest=earliest,
                                                          latest=latest)

        all_project_numbers = houston_permit_data['PROJECT_NO'].unique().tolist()

        # Step 2: Exclude any project numbers for which we already have data
        try:
            existing_data = pd.read_csv(target_path, index_col=0, header=None)
            existing_data.index = [str(id) for id in existing_data.index]
            all_project_numbers = [n for n in all_project_numbers if n not in existing_data.index.tolist()]
        except FileNotFoundError:  # In case we don't already have data
            print('Note that no pre-existing Houston permit approval data exists - starting from scratch.')
            pass

        # Step 3: Scrape new stuff
        print('Starting to scrape now')
        number_processes = 10
        result = speedy_get_project_statuses(all_project_numbers, number_processes=number_processes)
        result.to_csv(backup_path, mode='a')