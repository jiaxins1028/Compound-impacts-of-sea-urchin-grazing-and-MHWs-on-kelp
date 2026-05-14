###############################################################################
## This script extracts and processes E. radiata percent
## cover data from the Reef Life Survey (RLS) programme,
## then combines it with the existing ATRC time series to
## produce a single merged dataset used in downstream analyses.
##
##  Data sources:
##  ATRC  – pre-processed survey data available at:
##          github.com/jiaxins1028/MHW-impacts-on-Kelp-Populations/tree/main/Data
##          (loaded directly from er_atrc_id.csv; no extraction code needed here)
##  RLS   – https://catalogue-imos.aodn.org.au/geonetwork/srv/eng/catalog.search#/metadata/5a94a2cf-0810-44ea-a28d-cc8a5c30fbd7
##          enquiries@reeflifesurvey.com
###############################################################################

library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr) 
library(tidyr)
library(viridis)

#=======================RLS =============================== ---------------------------
# import data -------------------------------------------------------------
k_rls<-read.csv("...") #import RLS PQ seaweed data

MI.atrc<- alg.rls %>% filter(program =="RLS") #Filter just the RLS data

#convert characters to factors
k_rls$location<-as.factor(k_rls$location)
k_rls$country<-as.factor(k_rls$country)
k_rls$site_code<-as.factor(k_rls$site_code)
k_rls$RLS_category<-as.factor(k_rls$RLS_category)
k_rls$survey_id<-as.factor(k_rls$survey_id)
k_rls$dataset_id<-as.factor(k_rls$dataset_id)
k_rls$percent_cover<-as.integer(k_rls$percent_cover)
str(k_rls)

#convert survey_date from Character to date.
k_rls$survey_date <- as.Date(k_rls$survey_date, format = "%d/%m/%Y")
#Create new factor with survey year
k_rls$survey_year<-format(as.Date(k_rls$survey_date, format="%Y-%m-%d"),"%Y")
str(k_rls)

# Determine the unique time period for each location
location_time_periods_sites <- k_rls %>%
  group_by(location) %>%
  summarize(unique_category = list(unique(RLS_category)), unique_ids = list(unique(survey_id))) #unique_years = list(unique(survey_year)), unique_sites = list(unique(site_name)), 

complete_data <- location_time_periods_sites %>%
  rowwise() %>%
  mutate(data_frame = list(expand_grid(RLS_category = unique_category, survey_id = unique_ids))) %>%  # survey_year = unique_years, site_name = unique_sites, 
  select(-unique_category, -unique_ids) %>%  # -unique_years, -unique_sites, 
  unnest(cols = data_frame)

## Join the complete dataset with the original data
k_rls$survey_year<-as.factor(k_rls$survey_year)

k_rls_short = k_rls[, c("survey_id", "dataset_id", "location", "site_code", "site_name",  "survey_year",
                        "latitude", "longitude", "survey_date", "RLS_category", "percent_cover")]
data_joined <- complete_data %>%
  left_join(k_rls_short, by = c("location","survey_id", "RLS_category"))  #"survey_year", , "site_name"

df_joined_filled <- data_joined %>%
  group_by(survey_id) %>%
  mutate(across(c(site_code, site_name, latitude, longitude, survey_year), 
                ~ if_else(is.na(.), first(na.omit(.)), .))) %>%
  ungroup() %>%
  # Replace NA in percentage column with 0
  mutate(percent_cover = replace_na(percent_cover, 0))

## select specifc species
finaldata <- df_joined_filled %>% filter(RLS_category =="Ecklonia radiata")

# Average into site/location level ----------------------------------------
####average data, first by quadrat, then by survey_id, then by site and location

quad.id<-ddply(finaldata,c("survey_id","location","site_code","site_name", "latitude", "longitude",
                           "survey_year","RLS_category"),summarise, survey.mean = mean(percent_cover))
quad.site<-ddply(quad.id,c("location","site_name","survey_year",
                           "RLS_category", "latitude", "longitude"),summarise, survey.mean = mean(survey.mean))
#average data by location
quad.loca <-ddply(quad.site,c("location","survey_year","RLS_category"),summarise, survey.mean = mean(survey.mean))

# neat data into csv ------------------------------------------------------
write.table(quad.id, "er_rls_id.csv", sep = ",", row.names = FALSE, col.names = TRUE)


# plotting ----------------------------------------------------------------
quad.site$survey_year<-as.numeric(as.character(quad.site$survey_year))
quad.loca$survey_year<-as.numeric(as.character(quad.loca$survey_year))

df <- quad.site %>% filter(location =="Batemans")
df_loc <- quad.loca %>% filter(location =="Batemans")

ggplot() + 
  geom_line(data = df, aes(y = survey.mean, x = survey_year, group = site_name), alpha = 0.4) +  
  geom_point(data = df, aes(y = survey.mean, x = survey_year, group = site_name), alpha = 0.4) + 
  geom_line(data = df_loc, aes(y = survey.mean, x = survey_year, group = 1), color = 'black', linewidth = 1.5) + 
  scale_x_continuous(limits=c(2008,2019), breaks=seq(2008,2019,2))+
  labs(title = "RLS - Batemans",
       x = "Survey Year",
       y = "Kelp cover percent") +
  theme_format+
  theme(legend.position = "none")  # Hides the entire legend


# filling gaps in ATRC using RLS ------------------------------------------
df_rls <- read.csv("er_rls_id.csv") #import kelp data in RLS at survey_id level
df_atrc <- read.csv("er_atrc_id.csv") #import kelp data in ATRC at survey_id level (http://github.com/jiaxins1028/MHW-impacts-on-Kelp-Populations/tree/main/Data)

# Specify the two locations to check
target_locations <- c("Batemans", "Jervis Bay")

# Extract rows for specific locations where survey_year is missing in df2
missing_rows <- df_rls %>%
  filter(location %in% target_locations) %>%
  anti_join(df_atrc, by = c("location", "survey_year"))

# Combine the results
df_atrcrls <- bind_rows(df_atrc, missing_rows) %>% arrange(location, site_name, survey_year)  # Optional: Sort the result

write.table(df_atrcrls, file="er_atrcrls_id.csv", sep=',', row.names=F, col.names=T)
