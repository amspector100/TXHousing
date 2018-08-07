# TXHousing

Analysis of housing regulations in Texas. This is currently a bit messy and is in the process of being reorganized :)

## Setup & Steps to Reproduce

###

### Data

## Structure and Documentation

Almost all of the Python code lives inside the TXHousing python package, which is structured into three subpackages:
 "utilities", "data_processing", and "analysis". The utilities and data_processing packages are full of helper functions which
help streamline the code in the analysis package. For documentation of every function in the TXHousing package, see
 {} (to be inserted).

R code lives in the larger repository.

## Use

Because the TXHousing package is structured into subpackages which rely on each other, trying to run any module from
any of the subpackages will probably raise ImportErrors. Instead, you should import the TXHousing pacakge from outside
the