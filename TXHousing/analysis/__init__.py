"""Subpackage which does the analysis - this is where the magic happens"""
all = ['choropleth', 'mastermaps', 'BindColorMap', 'misc_calcs', 'parcel_graphs', 'permit_graphs', 'suburbs', 'zoning_graphs']

from . import choropleth
from . import mastermaps
from . import BindColorMap
from . import misc_calcs
from . import parcel_graphs
from . import permit_graphs
from . import suburbs
from . import zoning_graphs