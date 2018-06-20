import requests
import numpy as np
import pandas as pd
from lxml import html

def scrape_code_name(code):
    url = 'http://www.neighborhoodlink.com/zip/' + str(code)
    page = requests.get(url)
    tree = html.fromstring(page.content)
    print(tree)


scrape_code_name(10803)
