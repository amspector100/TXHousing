library(tidyverse)
library(maptools)
library(ggmap)
library(sf)

source('block-processing-functions.R')


#setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data/")
tx_block_path <- "ACS_2016_5YR_BG_48_TEXAS.gdb"
ny_block_path <- "ACS_2016_5YR_BG_36_NEW_YORK.gdb"
ca_block_path <- "ACS_2016_5YR_BG_06_CALIFORNIA.gdb"
il_block_path <- "ACS_2016_5YR_BG_17_ILLINOIS.gdb"
ma_block_path <- "ACS_2016_5YR_BG_25_MASSACHUSETTS.gdb"

old_cols <- sprintf('B25075e%d', 2:27)
new_cols <- c('0-10K', '10-15K', '15-20K', '20-25K', '25-30K', '30-35K',
              '35-40K', '40-50K', '50-60K', '60-70K', '70-80K', '80-90K',
              '90-100K', '100-125K', '125-150K', '150-175K', '175-200K',
              '200-250K', '250-300K')