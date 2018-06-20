# setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data")

# Source: 2016 5-year ACS tables DP-4

library(dplyr)
library(forcats)
library(ggplot2)
library(magrittr)
library(readr)
library(tibble)
library(tidyr)
library(RColorBrewer)

# Globals/column names
old_hsg_cols <- sprintf('HC01_VC%02d', 14:22)
old_cols <- c('GEO.display-label', old_hsg_cols)
new_hsg_cols <- c('OneUnitDetached', 'OneUnitAttached', 'TwoUnits', 'ThreeUnits', 'FiveUnits', 'TenUnits', 'TwentyUnits', 'MobileHome', 'Vehicle')
new_cols <- c('Area', new_hsg_cols)
stopifnot(length(new_cols) == length(old_cols))


# Initial read of data
path <- "Typology_Comparison_ACS_16_5YR_DP04/ACS_16_5YR_DP04.csv"
data <- read_csv(path) %>%
  select(one_of(old_cols)) %>%
  rename_at(vars(old_cols), ~ new_cols) %>%
  gather('HousingType', 'Units', new_hsg_cols)

# Refine data
data_coarse <- data %>%
  # Collapse types of housing
  mutate(HousingType=
           fct_collapse(HousingType,
                        '1'=c('OneUnitDetached', 'OneUnitAttached'),
                        '2\u20134' = c('TwoUnits','ThreeUnits'),
                        '5\u201319'= c('FiveUnits', 'TenUnits'),
                         #'2\u201319'=c('TwoUnits', 'ThreeUnits', 'FiveUnits', 'TenUnits'),
                        '20+'='TwentyUnits',
                        'Impermanent'=c('MobileHome', 'Vehicle')
           )
  ) %>%
  filter(HousingType != 'Impermanent') %>% 
  filter(Area != 'Los Angeles County, California') %>%
  # Rename areas
  #mutate(Area = ifelse(Area == 'Los Angeles County, California', 'LA County', Area)) %>%
  tidyr::separate(Area, into = c('Area'), sep=' city') %>% 
  group_by(HousingType, Area) %>%
  summarize(Units=sum(Units)/1000) %>%
  ungroup() %>%
  mutate(HousingType = factor(HousingType, 
                              levels = c('1', '2\u20134', '5\u201319', '20+'))) %>%
  mutate(Area = factor(Area,
                       levels = c("Austin",  "Dallas", "Houston", "Boston", "Chicago", 
                       "New York", "San Francisco", "Los Angeles")#, 
                       #"LA County")
                       )
         ) %>%
  group_by(Area) %>%
  mutate(Units = 100*Units/sum(Units)) %>%
  ungroup()

# Graph
bar_width = 0.9
colors <- brewer.pal(9, "Blues")[4:9]
plot_common <- 
  list(aes(x=HousingType, y=Units, fill = HousingType), 
       geom_col(position=position_dodge(bar_width)),
       labs(x='Units per building', y='Percent of Total Units'),
       facet_wrap(. ~ Area),
       scale_fill_manual(values=colors),
       scale_y_continuous(expand = c(0,5))
       )
       
gg <- ggplot(data_coarse) + plot_common

#setwd("C:/Users/aspector/Documents/R Projects/TX Housing")
#svg('Figures/comparison-housing-typology-v1.svg', w=10, h=8)
#print(gg)
#dev.off()



