import os

# Set working directory to the data directory.
file_directory = os.path.dirname(os.path.abspath(__file__))
parent_directory = os.path.split(file_directory)[0]
os.chdir(parent_directory)

import TXHousing.analysis.suburbs as suburbs
suburbs.texas_job_centers(num_polygons = 3000)