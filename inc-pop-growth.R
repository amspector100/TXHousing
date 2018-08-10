#! /usr/bin/RSCRIPT

library(bit64)
library(data.table)
library(geosphere)
library(ggplot2)
library(stringr)
library(tibble)

# Assemble 2000 statistics.

nas <- c('(X)', '*****')

demo_2000 <- fread('data/Census/2000-demographics.csv', na.strings=nas)
demo_2000_selected <- demo_2000[,.(GEOID00=Id2, Geography, TotalPop2000=`Number; Total population`, WhitePop=`Number; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - White alone`)]
demo_2000_selected[,WhitePct:=100*WhitePop/TotalPop2000]

# NB: No need to include "Geography" column twice!

econ_2000 <- fread('data/Census/2000-economics.csv', na.strings=nas)
econ_2000_selected <- econ_2000[,.(GEOID00=Id2, MedHHInc=`Number; INCOME IN 1999 - Households - Median household income (dollars)`)]

dt_2000 <- merge(demo_2000_selected, econ_2000_selected, on=GEOID00)

pop_2016 <- fread('data/Census/2016-demographics.csv')
pop_2016_selected <- pop_2016[,.(GEOID10=Id2, TotalPop2016=`Estimate; Total:`)]

housing_2016 <- fread('data/Census/2016-housing.csv')
housing_2016_selected <- housing_2016[,.(GEOID10=Id2, TotalHousing2016=`Estimate; HOUSING OCCUPANCY - Total housing units`)]

dt_2016 <- merge(pop_2016_selected, housing_2016_selected, on=GEOID10)

# Reallocate 2016 population to 2000 census tracts.
# 2000 tracts that were merged in 2010 get proportions of the 2016 population commensurate with their proportions of the 2010 population.

crosswalk <- fread('data/Census/tx-census-tract-crosswalk.csv')
dt_2016_with_crosswalk <- merge(crosswalk, dt_2016, by='GEOID10')
dt_2016_reduced <- dt_2016_with_crosswalk[,.(TotalPop2016=sum(TotalPop2016*POPPCT10/100), TotalHousing2016=sum(TotalHousing2016*HUPCT10/100)), by=GEOID00]

dt_demos <- merge(dt_2016_reduced, dt_2000, by='GEOID00')
gazetteer <- fread('data/Census/ustracts2k.csv')
gazetteer <- gazetteer[State == 'TX']
gazetteer[,LandAreaAcres:=LandAreaSqMeters/4046.856]
dt_final <- merge(dt_demos, gazetteer, by='GEOID00')

dt_final[,HousingGrowthPerAcre:=(TotalHousing2016-TotalHousing2000Gazetteer)/LandAreaAcres][,PopGrowthPerAcre:=(TotalPop2016-TotalPop2000Gazetteer)/LandAreaAcres]

# Split geo IDs into components
dt_final[,GeoIdStr:=str_pad(as.character(GEOID00), 11, 'left', '0')
][,StateFips:=as.integer(substr(GeoIdStr, 1, 2))
][,County:=as.integer(substr(GeoIdStr, 3, 5))
][,TractGroup:=as.integer(substr(GeoIdStr, 6, 9))
][,Tract:=ifelse(TractGroup %% 100 == 0,
as.character(TractGroup),
paste(TractGroup, substr(GeoIdStr, 10, 11), sep='.'))]

# Modifications from the initial 3/28 script begin here.
library(sf)
library(magrittr)
library(dplyr)

# Start by getting geodata for census tracts in 2000
tracts_2000 <- sf::st_read("data/Census/Texas_Tract_Boundaries_2000/tl_2010_48_tract00.shp") %>%
  dplyr::mutate(CTIDFP00 = as.numeric(as.character(CTIDFP00))) %>%
  dplyr::rename(GEOID00 = CTIDFP00)
dt_final <- as.tibble(dt_final) %>%
  dplyr::mutate(GEOID00 = as.numeric(GEOID00))
geodata_2000 <- inner_join(tracts_2000, dt_final, by = c('GEOID00' = 'GEOID00'))

metro_area_names <- c('Austin-Round Rock, TX', 
                      'Dallas-Fort Worth-Arlington, TX',
                      'Houston-The Woodlands-Sugar Land, TX')
place_names <- c('Austin', 'Dallas', 'Houston')

# Now get geodata for metro-areas (CBSA) and for places 
metro_shapes <- st_read("data/Census/cb_2017_us_cbsa_500k/cb_2017_us_cbsa_500k.shp") %>%
  dplyr::mutate(NAME = as.character(NAME)) %>%
  dplyr::filter(NAME %in% metro_area_names)
metro_shapes <- st_cast(metro_shapes, 'POLYGON')
places <- sf::st_read("data/Census/cb_2017_48_place_500k/cb_2017_48_place_500k.shp") %>%
  dplyr::mutate(NAME = as.character(NAME)) %>%
  dplyr::filter(NAME %in% place_names)
places <- st_cast(places, "POLYGON")

# Subset to metro areas/places
# These intersections are a bit slow because this is R
# and we're not using spatial trees.
final_columns <- c('GEOID00', 'COUNTYFP00', 'Tract', 'MedHHInc',
                   'PopGrowthPerAcre', 'WhitePct', 'NAME', 'Geography')
core_areas <- st_intersection(geodata_2000, places) %>%
  mutate('Geography' = 'Core City') %>%
  select(final_columns)
metro_areas <- st_intersection(geodata_2000, metro_shapes) %>%
  mutate('Geography' = 'Metropolitan Area') %>%
  select(final_columns)

# Alternatively...
final_data <- as.data.frame(rbind(core_areas, metro_areas)) %>%
  mutate(City = case_when(
    grepl('Austin', NAME) ~ 'Austin',
    grepl('Dallas', NAME) ~ 'Dallas',
    grepl('Houston', NAME) ~ 'Houston',
    TRUE ~ 'Unknown'
  ))

# Plot number 1: White demographic
gg <- ggplot(final_data, aes(x=WhitePct, y=PopGrowthPerAcre, color=City)) +
  geom_point() + 
  labs(title = 'Population Growth versus Percentage of White Individuals',
       xlab = 'Percentange of White Individuals (2000)',
       ylab = 'Population Growth per Acre in Census Tracts (2000-2016)',
       caption = 'Data from 2000 Census and 2016 ACS') +
  facet_grid(Geography ~ City)

svg('Figures/whitepct_v_popgrowth.svg', width=10, height=8)
print(gg)
dev.off()

# Plot number 2: Income demographic
gg2 <- ggplot(final_data, aes(x=MedHHInc, y=PopGrowthPerAcre, color=City)) +
  geom_point() + 
  labs(title = 'Population Growth versus Median Household Income',
       xlab = 'Median Household Income (2000)',
       ylab = 'Population Growth per Acre in Census Tracts (2000-2016)',
       caption = 'Data from 2000 Census and 2016 ACS') +
  facet_grid(Geography ~ City) 

  
svg('Figures/income_v_popgrowth.svg', width=10, height=8)
print(gg2)
dev.off()
