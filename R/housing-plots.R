#! /usr/bin/RSCRIPT

# Sources: 2000 SF3 and 2010, 2016 5-year ACS tables DP-4

library(dplyr)
library(forcats)
library(ggplot2)
library(magrittr)
library(readr)
library(tibble)
library(tidyr)

# Columns to preserve and their new names.
old_hsg_cols_2000 <- sprintf('HC01_VC%02d', 3:11)
old_cols_2000 <- c('GEO.display-label', old_hsg_cols_2000)

old_hsg_cols_2010 <- sprintf('HC01_VC%02d', 13:21)
old_cols_2010 <- c('GEO.display-label', old_hsg_cols_2010)

old_hsg_cols_2016 <- sprintf('HC01_VC%02d', 14:22)
old_cols_2016 <- c('GEO.display-label', old_hsg_cols_2016)

new_hsg_cols <- c('OneUnitDetached', 'OneUnitAttached', 'TwoUnits', 'ThreeUnits', 'FiveUnits', 'TenUnits', 'TwentyUnits', 'MobileHome', 'Vehicle')
new_cols <- c('Area', new_hsg_cols)
stopifnot(length(new_cols) == length(old_cols_2000))
stopifnot(length(new_cols) == length(old_cols_2010))
stopifnot(length(new_cols) == length(old_cols_2016))

# Read in and process data. Filtering precedes row-binding to save memory.
data_2000 <- read_csv('data/Census/DEC_00_SF3_DP4.csv') %>%
  select(one_of(old_cols_2000)) %>%
  rename_at(vars(old_cols_2000), ~ new_cols) %>%
  mutate(Year='2000 (Decennial Census)')

data_2010 <- read_csv('data/Census/ACS_10_5YR_DP04.csv') %>%
  select(one_of(old_cols_2010)) %>%
  rename_at(vars(old_cols_2010), ~ new_cols) %>%
  mutate(Year='2010 (five-year ACS)')

data_2016 <- read_csv('data/Census/ACS_16_5YR_DP04.csv') %>%
  select(one_of(old_cols_2016)) %>%
  rename_at(vars(old_cols_2016), ~ new_cols) %>%
  mutate(Year='2016 (five-year ACS)')

years <- c('2000 (Decennial Census)', 
           '2010 (five-year ACS)', 
           '2016 (five-year ACS)')

data <- bind_rows(data_2000, data_2010, data_2016) %>%
  gather(HousingType, Units, -Area, -Year) %>%
  mutate(Units=as.integer(Units),
         HousingType=parse_factor(HousingType, new_hsg_cols),
         Year=
           parse_factor(Year, years))

data_coarse <- data %>%
  # Collapse types of housing
  mutate(HousingType=
         fct_collapse(HousingType,
                      '1'=c('OneUnitDetached', 'OneUnitAttached'),
                      # '2\u20134'=c('TwoUnits', 'ThreeUnits'),
                      # '5\u20139'='FiveUnits',
                      # '10\u201319'='TenUnits',
                      '2\u201319'=c('TwoUnits', 'ThreeUnits', 'FiveUnits', 'TenUnits'),
                      '20+'='TwentyUnits',
                      'Impermanent'=c('MobileHome', 'Vehicle')
                      )
        ) %>%
  filter(HousingType != 'Impermanent') %>% 
  group_by(HousingType, Year, Area) %>%
  summarize(Units=sum(Units)/1000) %>%
  ungroup()

# Compute percent changes and create bar chart labels.
# Labels are 1000s of units in 2000, 1000s of units + % change in 2015.
label_data <- function(df) {
  df2 <- df %>%
    arrange(Year) %>%
    ungroup() %>%
    group_by(HousingType, Area) %>%
    mutate(PercentChange=100*(Units-lag(Units))/lag(Units),
           Label=ifelse(Year == '2000 (Decennial Census)',
                        sprintf("%.0f", Units),
                        sprintf("%.0f (%.1f%%)", Units, PercentChange))) %>%
    select(-PercentChange) %>%
    ungroup()
}
# Make a 3x3 plot of housing trends.

# First, compile core city data.
core_cities <- c('Austin', 'Dallas', 'Houston')
core_cities_census <- sprintf('%s city, Texas', core_cities)
core_cities_data <- data_coarse %>%
  filter(Area %in% core_cities_census) %>%
  mutate(Area=factor(Area, levels=core_cities_census, labels=core_cities),
         Scale='Core city') %>%
  label_data()

# Next, a list of all counties in Texas.
regions <- unique(data_coarse$Area)
tx_counties_census <- regions[grepl("County, Texas", regions, fixed=T)]
tx_counties <- gsub(" County, Texas", "", tx_counties_census)

# Shorten county names for county-level data.
counties_data <- data_coarse %>%
  filter(Area %in% tx_counties_census) %>%
  mutate(Area=factor(Area, levels=tx_counties_census, labels=tx_counties))

core_counties_data <- counties_data %>%
  filter(Area %in% c('Dallas', 'Travis', 'Harris')) %>%
  # Fun fact: the magrittr pipe works on any object, not just data frames.
  mutate(Area=fct_drop(Area) %>%
              fct_relevel(c('Travis', 'Dallas', 'Harris')) %>%
              fct_recode('Austin'='Travis','Houston'='Harris'),
         Scale='Core county') %>%
  label_data()

# Last: metropolitan areas, defined by counties
metro_areas_data <- counties_data %>%
  mutate(Area=
         fct_collapse(Area,
                      'Austin2'=c('Travis', 'Williamson'),
                      # There is a distinct Austin County! Using 'Austin' as a separate factor name would get them confused.
                      'Dallas'=c('Dallas', 'Tarrant', 'Denton',
                              'Collin', 'Rockwall'),
                      'Houston'=c('Harris', 'Montgomery', 'Galveston',
                                'Fort Bend', 'Brazoria'),
                      )
         ) %>%
  filter(Area %in% c('Austin2', 'Dallas', 'Houston')) %>%
  mutate(Area=fct_recode(Area, Austin='Austin2') %>%
              fct_drop() %>%
              fct_relevel('Austin', 'Dallas', 'Houston')) %>%
  group_by(Year, Area, HousingType) %>%
  summarize(Units=sum(Units)) %>%
  ungroup() %>%
  label_data() %>%
  mutate(Scale='Metropolitan area')

combined_data <-
  bind_rows(core_cities_data, core_counties_data, metro_areas_data) %>%
  mutate(Scale=
         factor(Scale,
                levels=c('Core city', 'Core county', 'Metropolitan area')))


# ggplot theme elements
bar_width <- 0.95
plot_common <-
  list(aes(x=HousingType, y=Units, fill=Year, group=Year, label=Label),
       geom_col(position=position_dodge(bar_width)),
       geom_text(position=position_dodge(bar_width),
                 size = 3, angle = 90, hjust=-0.055, vjust = -0.2),
       scale_fill_manual(values=c('#45505f', '#4ab1e4', '#000000')),
       theme_bw(),
       labs(x='Units per building', y='Total units (thousands)'),
       theme(legend.position=c(1, 1),
             legend.justification=c(1, 1),
             legend.box.margin=margin(10, 10, 10, 10),
             legend.title=element_blank(),
             axis.ticks.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y=element_line(color='gray80'),
             strip.text=element_text(color='black'),
             axis.text=element_text(color='black', size = 15),
             strip.background=element_blank()),
       facet_grid(Scale ~ Area),
       # Asymmetrical top and bottom margins
       scale_y_continuous(expand=expand_scale(mult=c(0, .2)))
       )

gg <- ggplot(combined_data) + plot_common
svg('Figures/texas-housing-typology.svg', height=10, width=10)
print(gg)
dev.off()
