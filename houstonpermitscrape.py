from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup
import time
import pandas as pd
from multiprocessing import Pool
import sys

# The file to save everything to
target_path = 'shared_data/houston_permit_statuses.csv'
backup_path = 'shared_data/houston_permit_statuses_backup.csv'


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
            print('This process has finished {} ids'.format(counter))

        # Some very crude error handling
        try:
            # Go to website
            browser.get(url)

            # Search for specific project
            input_bar = browser.find_element_by_name('ProjectNo')
            input_bar.send_keys(project_number)
            input_bar.send_keys(Keys.RETURN) # Hit return after sending text

            # Get content and return result
            content = browser.find_elements_by_class_name('content')
            content_list = [item.text for item in content]
        except:
            print("Unexpected error for {}".format(project_number), sys.exc_info()[0])
            return result_dictionary

        if 'These plans have NOT been approved for permitting' in content_list:
            result_dictionary[project_number] = 0.0 # Not approved
        elif len(content_list) <= 3:
            result_dictionary[project_number] = float('NaN') # No data
        else:
            result_dictionary[project_number] = 1.0 # Approved

        counter += 1

    # Quit browser
    browser.quit()
    result_dictionary = pd.Series(result_dictionary)
    return result_dictionary

# Use multithreading to do the same thing
def speedy_get_project_statuses(all_project_numbers, number_processes = 3):

    time0 = time.time()
    chunksize = round(len(all_project_numbers)/number_processes) + 1

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
    results.to_csv(target_path, mode = 'a', header = False)
    print('Took {} seconds for {} project numbers'.format(time.time() - time0, len(all_project_numbers)))
    return results


if __name__ == '__main__':

    from helpers import process_houston_permit_data

    # Step 1: Get all the unique project numbers
    #searchfor = ['NEW S.F.', 'NEW SF', 'NEW SINGLE', 'NEW TOWNHOUSE', 'NEW AP', 'NEW HI-']
    searchfor = [''] # Almost everything is residential demolition, and anyway, commercial demos are included in Dallas
    searchin = ['PROJ_DESC']
    #kind = 'structural'
    kind = 'demolition'
    earliest = 2010
    latest = None

    print('Processing permit data')
    houston_permit_data = process_houston_permit_data(searchfor = searchfor,
                                                      searchin = searchin,
                                                      kind = kind,
                                                      earliest = earliest,
                                                      latest = latest)

    all_project_numbers = houston_permit_data['PROJECT_NO'].unique().tolist()

    # Step 2: Exclude any project numbers for which we already have data
    try:
        existing_data = pd.read_csv(target_path, index_col = 0, header = None)
        existing_data.index = [str(id) for id in existing_data.index]
        all_project_numbers = [n for n in all_project_numbers if n not in existing_data.index.tolist()]
    except FileNotFoundError: # In case we don't already have data
        print('Note that no pre-existing data exists')
        pass

    # Step 3: Scrape new stuff
    print('Starting to scrape now')
    number_processes = 10
    result = speedy_get_project_statuses(all_project_numbers, number_processes = number_processes)
    result.to_csv(backup_path, mode = 'a')