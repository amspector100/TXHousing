import os

# Set working directory to the data directory.
file_directory = os.path.dirname(os.path.abspath(__file__))
parent_directory = os.path.split(file_directory)[0]
os.chdir(parent_directory)

import TXHousing.data_processing.zoning as zoning
import TXHousing.analysis.suburbs as suburbs
suburbs.analyze_transportation_networks(names = ['Austin', 'Dallas', 'Houston'],
                                        zoning_inputs = [zoning.austin_inputs, zoning.dallas_inputs, zoning.houston_inputs])
