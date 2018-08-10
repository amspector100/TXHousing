# Wrapper which generates all of the graphs which were made in R

# This one is fast (~10 seconds)
source('R/inc-pop-growth.R')
rm(list = ls())

# Also fast (~5 sec)
source('R/housing-plots.R')
rm(list = ls())

# Also fast (~5 sec)
source('R/comparison-housing-plots.R')
rm(list = ls())

# This one takes 2-3 minutes or so
source('R/density-histogram.R')
rm(list = ls())