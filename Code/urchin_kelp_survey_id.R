library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr) 
library(tidyr)
library(viridis)
library(arsenal)


# import data -------------------------------------------------------------
### Raw NRMN seaweed data can be downloaded from Raw data and instructions can be accessed from 
### https://catalogue-imos.aodn.org.au/geonetwork/srv/eng/catalog.search#/metadata/ec424e4f-0f55-41a5-a3f2-726bc4541947
### Read https://github.com/jiaxins1028/MHW-impacts-on-Kelp-Populations/blob/main/rawdata_download_use.md

alg.rls<-read.csv("...") #import raw urchin data
alg.rls <- mutate(alg.rls, number=round(total/50,3))  

MI.atrc<- alg.rls %>% filter(program =="ATRC") #Filter ATRC / RLS data
#convert characters to factors
MI.atrc$location<-as.factor(MI.atrc$location)
MI.atrc$area<-as.factor(MI.atrc$area)
MI.atrc$country<-as.factor(MI.atrc$country)
MI.atrc$ecoregion<-as.factor(MI.atrc$ecoregion)
MI.atrc$realm<-as.factor(MI.atrc$realm)
MI.atrc$site_code<-as.factor(MI.atrc$site_code)
MI.atrc$phylum<-as.factor(MI.atrc$phylum)
MI.atrc$class<-as.factor(MI.atrc$class)
MI.atrc$order<-as.factor(MI.atrc$order)
MI.atrc$family<-as.factor(MI.atrc$family)

str(MI.atrc)

#convert survey_date from Character to date.
MI.atrc$survey_date <- as.Date(MI.atrc$survey_date, format = "%Y-%m-%d")
#Create new factor with survey year
MI.atrc$survey_year<-format(as.Date(MI.atrc$survey_date, format="%d-%m-%YY"),"%Y")
# MI.atrc$survey_year<-as.factor(MI.atrc$survey_year)
str(MI.atrc)


## Determine the unique time period for each location
location_time_periods_sites <- MI.atrc %>%
  group_by(location) %>%
  summarize(unique_taxon = list(unique(taxon)), unique_ids = list(unique(survey_id))) #unique_years = list(unique(survey_year)), unique_sites = list(unique(site_name)), 

complete_data <- location_time_periods_sites %>%
  rowwise() %>%
  mutate(data_frame = list(expand_grid(taxon = unique_taxon, survey_id = unique_ids))) %>%  # survey_year = unique_years, site_name = unique_sites, 
  select(-unique_taxon, -unique_ids) %>%  # -unique_years, -unique_sites, 
  unnest(cols = data_frame) 

## Join the complete dataset with the original data
MI.atrc$survey_year<-as.factor(MI.atrc$survey_year)
MI.atrc_short = MI.atrc[, c("survey_id", "location", "site_code", "site_name",
                            "latitude", "longitude", "survey_date", "program", 
                            "taxon", "number", "survey_year", "size_class", 'block')]

data_joined <- complete_data %>%
  left_join(MI.atrc_short, by = c("location","survey_id", "taxon"))  #"survey_year", , "site_name"


df_joined_filled <- data_joined %>%
  group_by(survey_id) %>%
  mutate(across(c(site_code, site_name, latitude, longitude, program, survey_year, survey_date), ~ if_else(is.na(.), first(na.omit(.)), .))) %>%
  ungroup() %>%
  # Replace NA in number column with 0
  mutate(number = replace_na(number, 0))

finaldata <- df_joined_filled %>% filter(taxon == "Centrostephanus rodgersii")   


# Average into site/location level ----------------------------------------
####average data, first by survey ID, then by site/location
quad.size<-ddply(finaldata,c("survey_id","location","site_name","survey_year", "latitude", "longitude", 
                             "survey_date", "taxon", 'block'),summarise, number = sum(number))
quad.block<-ddply(quad.size,c("survey_id","location","site_name","survey_year", "latitude", "longitude", 
                              "survey_date", "taxon"),summarise, number = mean(number))
quad.id<-ddply(quad.block,c("survey_id","location","site_name","survey_year", "latitude", "longitude", 
                            "survey_date"),summarise, number = sum(number))

quad.site<-ddply(quad.id,c("location","site_name","survey_year"),summarise, number = mean(number))
#average data by location
quad.loc<-ddply(quad.site,c("location","survey_year"),summarise, number = mean(number))
# str(quad.site1)


# neat data into csv ------------------------------------------------------
write.table(quad.id, file="cr_atrc_id.csv", sep=',', row.names=F, col.names=T)  ## or RLS data -- cr_rls_id.csv



# filling gaps in ATRC using RLS ------------------------------------------
df_rls <- read.csv("cr_rls_id.csv") #import urchin data in RLS at survey_id level
df_atrc <- read.csv("cr_atrc_id.csv") #import urchin data in RLS at survey_id level

# Specify the two locations to check
target_locations <- c("Batemans", "Jervis Bay")

# Extract rows for specific locations where survey_year is missing in df2
missing_rows <- df_rls %>%
  filter(location %in% target_locations) %>%
  anti_join(df_atrc, by = c("location", "survey_year"))

# Combine the results
df_atrcrls <- bind_rows(df_atrc, missing_rows) %>% arrange(location, site_name, survey_year)  # Optional: Sort the result

write.table(df_atrcrls, file="cr_atrcrls_id.csv", sep=',', row.names=F, col.names=T)



# combine comparable kelp and urchin survey_ids -------------------------------------------------------------

df_kelp<-read.csv("er_atrcrls_id.csv") #import kelp data (ATRC + RLS; Output from `kelp_atrc_rls_surveyid.R`)
df_urchin<-read.csv("cr_atrcrls_id.csv") #import kelp data (ATRC + RLS)


# merge kelp and urchin dfs
combined_df <- df_kelp %>%
  inner_join(df_urchin, by = "survey_id")

# find duplicated survey_ids
combined_df$survey_id[duplicated(combined_df$survey_id)]

# neat data into csv ------------------------------------------------------
write.table(combined_df, "ercr_atrcrls_id.csv", sep = ",", row.names = FALSE, col.names = TRUE)


