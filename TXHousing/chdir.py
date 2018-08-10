import os

# All this file does is set working directory to the data directory. It's only used in the testing modules.
file_directory = os.path.dirname(os.path.abspath(__file__))
parent_directory = os.path.split(file_directory)[0]
os.chdir(parent_directory)
