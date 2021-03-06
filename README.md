# TXHousing

Analysis of housing regulations in Texas, hosted on this repo in the aim of making the research reproducible. 
This repo is mostly finished but still under a bit of development. If you have questions, check out the 
documentation [linked here](https://amspector100.github.io/TXHousing) or email other questions to amspector100@gmail.com.

## Structure and Documentation

All of the Python code lives inside the TXHousing python package, which is structured into three subpackages:
"utilities", "data_processing", and "analysis". The utilities and data_processing packages are full of helper functions
which help streamline the code in the analysis package. For documentation of the TXHousing package, see
[amspector100.github.io/TXHousing](https://amspector100.github.io/TXHousing). 

We did use some R code, which lives in the R subdirectory. It's very simple, so it's not packaged and there's no official 
documentation site - but there might be in the future! (It might also eventually migrate to Python).

## General Requirements

This project basically draws upon most of the geospatial stack in R/Python. For a more specific list of requirements,
see the requirements.txt file within the TXHousing subdirectory and the (relatively primitive) DESCRIPTION file in the
R subdirectory.

## Setup & Steps to Reproduce

### Clone

Clone the repo from GitHub, i.e.:

    git clone https://github.com/amspector100/TXHousing

### Data

The data for the project is divided into a 'data' directory, which contains almost all of the raw data and the processed 
parcel data, and a 'shared_data' directory, which contains some smaller caches, calculations, and a couple of  raw data 
files which are not available on the internet. **Unfortunately, the project required too much raw data to store on 
GitHub, so only the shared_data directory is available directly in this repo.** However, you can download the data from
 the following two links on google drive: [data](https://drive.google.com/open?id=1ergEpaS-vf7R-Pt_GX0yxVpYvjPxkVTF), 
 [shared_data](https://drive.google.com/open?id=1P7sBfr1suDWvi6k7LooFwDTGajvC8U1I).

**Alternatively**, if you want to reproduce everything from scratch, the testing modules will alert you to what data you are missing
at any given time and also provide you a link at which you may download the data. Basically, you should run the 
data_processing testing module (as described later) and download data until the data_processing module stops throwing 
errors. After this, you should have all of the raw data of the project with the exception of parcel files used to 
do analysis of the suburbs. The paths and errors for these are listed separately in the parcel.py module in the 
data_processing subpackage, the reason for this being that the suburbs data is enormous and ugly, and they're not really 
required to do most of the analysis in the project. Note that once you download and extract the data, you will have to
adjust the paths to match the global variables supplied in the data_processing subpackage. 

All data which is not available on the internet is stored in the 'shared_data' directory for this reason.

### Getting started

After cloning the package and setting up the data, you should first run the two test modules, test_processing and 
test_utilities. These are **not** comprehensive testing suites but will let you know if anything is going
catastrophically wrong inside the utilities or data processing subpackages. To run them, you should simply navigate to the 
test directory inside the repo and run the files, using a command like:

    python test_utilities.py
    python test_processing.py

If you have not downloaded the data directory, you will need to generate a couple of data caches before running analysis.
You can generate them by running the R/tarrant-parcel-processing.R script and the generate_caches.py script, **in that order.**
Note that generating the caches for parcel data and houston permit statuses will likely take several hours and a lot of 
memory. None of this is necessary if you have downloaded the [data](https://drive.google.com/open?id=1ergEpaS-vf7R-Pt_GX0yxVpYvjPxkVTF) directory.

If you want to reproduce all the graphs, you should be able to run the 'generate-all-graphs.R' and the 'generate_graphs.py'
scripts after generating the caches. You can also call functions in the TXHousing.analysis.misc_calcs module to reproduce
any other miscellaneous calculations. 

### Notes on More Detailed Use

Because the TXHousing package is structured into subpackages which rely on each other, trying to run any module from
any of the subpackages will probably raise ImportErrors. Instead, you should import the TXHousing package either from 
the repository or from a folder above the repository. 

For very specific documentation of the TXHousing package, see
[amspector100.github.io/TXHousing](https://amspector100.github.io/TXHousing). 