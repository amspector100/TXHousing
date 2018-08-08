# TXHousing

Analysis of housing regulations in Texas. This repo, and this README, are still under construction :)

## Setup & Steps to Reproduce

### Clone

### Data

## Structure and Documentation

All of the Python code lives inside the TXHousing python package, which is structured into three subpackages:
"utilities", "data_processing", and "analysis". The utilities and data_processing packages are full of helper functions
which help streamline the code in the analysis package. For documentation of the TXHousing package, see
[amspector100.github.io/TXHousing](amspector100.github.io/TXHousing).

R code lives in the larger repository.

## Use

Because the TXHousing package is structured into subpackages which rely on each other, trying to run any module from
any of the subpackages will probably raise ImportErrors. Instead, you should import the TXHousing pacakge from outside
the