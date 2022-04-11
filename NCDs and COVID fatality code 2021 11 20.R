# This is the analytical code for the manuscript entitled 'Burden of Non-Communicable Diseases (NCDs) and COVID-19 Fatality Ratio: An Ecological Study in 104 Low- and Middle-Income Countries (LMICs' 
# The analyses relies on the following datasets: 

# Our World in Data on COVID-19 - available from: https://covid.ourworldindata.org/data/owid-covid-data.csv
# The World Bank - available from: https://data.worldbank.org/
# The Global Burden of Disease are not publicly available; full data for this analysis are not therefore provided.

# The code was prepared by Peter Tennant using R 4.1.0  

#### PREPARATION ####

#Load and/or install required packages

list.of.packages <- c("MASS", "tidyverse", "dplyr", "readr", "ggplot2", "forecast", "zoo", "bestglm", "rJava", "glmulti")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(MASS); library(tidyverse); library(dplyr); library(readr); library(ggplot2); library(forecast); library(zoo); library(bestglm); library(rJava); library(glmulti)

rm(list.of.packages, new.packages)

#Define list of 113 LMICs for examination
List_countryid  <- list("AFG", "ALB", "DZA", "AGO", "ARG", "ARM", "AZE", "BGD", "BLR", "BLZ", "BEN", "BTN", "BOL", "BIH", "BWA", "BRA", "BGR", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "CHN", "COL", "COM", "COD", "COG", "CRI", "CIV", "CUB", "DJI", "DOM", "ECU", "EGY", "SLV", "GNQ", "ERI", "SWZ", "ETH", "FJI", "GAB", "GMB", "GEO", "GHA", "GTM", "GIN", "GNB", "GUY", "HTI", "HND", "IND", "IDN", "IRN", "IRQ", "JAM", "JOR", "KAZ", "KEN", "KGZ", "LBN", "LSO", "LBR", "LBY", "MDG", "MWI", "MYS", "MDV", "MLI", "MRT", "MEX", "MNG", "MNE", "MAR", "MOZ", "MMR", "NAM", "NPL", "NIC", "NER", "NGA", "MKD", "PAK", "PNG", "PRY", "PER", "PHL", "RUS", "RWA", "STP", "SEN", "SRB", "SLE", "SOM", "ZAF", "LKA", "SDN", "SUR", "SYR", "TJK", "TZA", "THA", "TGO", "TUN", "TUR", "UGA", "UKR", "UZB", "VEN", "VNM", "YEM", "ZMB", "ZWE")

#Import 'outcome' (COVID) data and reduce to relevant list of countries
COVID_data      <- read_csv("https://github.com/owid/covid-19-data/raw/master/public/data/owid-covid-data.csv", show_col_types = FALSE)
COVID_data      <- rename(COVID_data, countryid = iso_code)
COVID_data      <- subset(COVID_data, countryid %in% List_countryid)

#Reduce to relevant variables:
COVID_data <- COVID_data %>% 
  select(c(countryid, continent, location, date, total_cases, new_cases, total_deaths, new_deaths, reproduction_rate, new_tests, total_tests, new_vaccinations, total_vaccinations, stringency_index, population, median_age, aged_65_older, aged_70_older, gdp_per_capita, hospital_beds_per_thousand, life_expectancy, human_development_index)) %>% 
  group_by(countryid) %>%
  mutate(day_in_data      = row_number()) %>% 
  ungroup(countryid)

# Interpolate missing values:
# Create a blank dataframe to save smoothed values
COVID_smoothed <- tibble(data.frame(countryid = character(0), day_in_data=numeric(0),  s_stringency_index=numeric(0), stringsAsFactors=F))

#Create list of countries requiring interpolation for stringency and reproduction rate:
List_imputing  <- list("AFG", "ALB", "DZA", "AGO", "ARG", "AZE", "BGD", "BLR", "BLZ", "BEN", "BTN", "BOL", "BIH", "BWA", "BRA", "BGR", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "CHN", "COL", "COD", "COG", "CRI", "CIV", "CUB", "DJI", "DOM", "ECU", "EGY", "SLV", "ERI", "SWZ", "ETH", "GAB", "GMB", "GEO", "GHA", "GTM", "GIN", "GUY", "HTI", "HND", "IND", "IDN", "IRN", "IRQ", "JAM", "JOR", "KAZ", "KEN", "KGZ", "LBN", "LSO", "LBR", "LBY", "MDG", "MWI", "MYS", "MLI", "MRT", "MEX", "MNG", "MAR", "MOZ", "MMR", "NAM", "NPL", "NIC", "NER", "NGA", "PAK", "PNG", "PRY", "PER", "PHL", "RUS", "RWA", "SEN", "SRB", "SLE", "SOM", "ZAF", "LKA", "SDN", "SUR", "SYR", "TJK", "TZA", "THA", "TGO", "TUN", "TUR", "UGA", "UKR", "UZB", "VEN", "VNM", "YEM", "ZMB", "ZWE")

lapply(List_imputing, function(COUNTRY) {
  
  Data_Time_Series <- ts(data = subset(COVID_data, countryid == COUNTRY))
  
  #Create smoothed values:
  imputed_stringency    <- na.interp(Data_Time_Series[,"stringency_index"], lambda = "auto")
  imputed_rr            <- na.interp(Data_Time_Series[,"reproduction_rate"], lambda = "auto")
  
  imputed_data          <- do.call(rbind, Map(data.frame, countryid=COUNTRY, day_in_data=Data_Time_Series[,"day_in_data"], s_stringency_index=imputed_stringency, s_reproduction_rate=imputed_rr))
  
  COVID_smoothed        <<- tibble(rbind(COVID_smoothed, imputed_data))
  
})

# Merge interpolated values into main dataset
COVID_smoothed        <- tibble(remove_rownames(COVID_smoothed))
COVID_data            <- left_join(COVID_data, COVID_smoothed, by=c("countryid", "day_in_data"))

#Remove imputation material:
rm(COVID_smoothed, List_imputing)

#Replace remaining missing values with zeros:
COVID_data$new_cases                <- ifelse(is.na(COVID_data$new_cases), 0, COVID_data$new_cases)
COVID_data$new_deaths               <- ifelse(is.na(COVID_data$new_deaths), 0, COVID_data$new_deaths)
COVID_data$new_deaths               <- ifelse(COVID_data$new_deaths<0, 0, COVID_data$new_deaths)
COVID_data$new_cases_per_million    <- (COVID_data$new_cases/COVID_data$population)*1000000

#Calculate measures of stringency and timeliness

COVID_data        <- COVID_data %>% 
  group_by(countryid) %>% 
  mutate(date = as.Date(date, format="%d/%m/%Y", origin="1970-01-01")) %>% 
  mutate(seven_day_average     = rollmean(new_cases_per_million*10, 7, fill="extend"),
         first_10_cases_flag   = if_else(total_cases>=10 & lag(total_cases<10, 1) & date<as.Date("2021-04-01"), 1, 0),
         first_10_cases_flag   = if_else(is.na(first_10_cases_flag), 0, first_10_cases_flag),
         first_lockdown_flag   = if_else( (s_stringency_index>=50 & lag(s_stringency_index<50)), 1, 0),
         first_lockdown_flag   = if_else( (s_stringency_index>=50 & day_in_data==1), 1, first_lockdown_flag),
         running_flag          = sum(first_lockdown_flag),
         days_at_10            = sum(if_else(first_10_cases_flag==1, day_in_data, 0)),
         days_at_lockdown      = sum(if_else(first_lockdown_flag==1 & running_flag==1, day_in_data, 0)),
         days_at_lockdown      = if_else(days_at_lockdown==0, max(day_in_data), days_at_lockdown),
         rate_at_lockdown      = round(sum(if_else(first_lockdown_flag==1, seven_day_average, 0)), digits=3),
         rate_at_lockdown      = if_else(rate_at_lockdown<0, 0, rate_at_lockdown),
         time_to_lockdown      = days_at_lockdown - days_at_10)

#restrict time period to 1st April 2021
COVID_data        <- COVID_data %>% 
  filter(date <= as.Date("2021-04-01")) 

#And for sensitivity analyses
COVID_data_Jan        <- COVID_data %>% 
  filter(date <= as.Date("2021-01-01")) 

COVID_data_Aug        <- COVID_data %>% 
  filter(date <= as.Date("2021-08-01")) 

#Summarise key data and collapse dataset
#Primary dataset

COVID_data <- COVID_data %>% 
  group_by(countryid) %>%
  mutate(lag_cases             = lag(total_cases, n=16),
         lag_cases             = ifelse(is.na(lag_cases)==TRUE, 0, lag_cases)) %>% 
  summarise(continent          = first(continent),
            location           = first(location),
            cases              = max(total_cases, na.rm=TRUE),
            lag_cases          = max(lag_cases, na.rm=TRUE),
            deaths             = max(total_deaths, na.rm=TRUE),
            tests              = max(total_tests, na.rm=TRUE),
            vaccinations       = max(total_vaccinations, na.rm=TRUE),
            stringency         = mean(s_stringency_index, na.rm=TRUE),
            reproduction       = mean(s_reproduction_rate, na.rm=TRUE),
            rate_at_lockdown   = first(rate_at_lockdown),
            time_to_lockdown   = first(time_to_lockdown),
            population         = first(population),
            gdp_per_capita     = first(gdp_per_capita),
            hospital_beds      = first(hospital_beds_per_thousand), 
            hdi                = first(human_development_index)) %>% 
  mutate(tests                 = ifelse(is.infinite(tests)==TRUE, NA, tests),
         vaccinations          = ifelse(is.infinite(vaccinations)==TRUE, NA, vaccinations))

#Sensitivity analysis

COVID_data_Jan <- COVID_data_Jan %>% 
  group_by(countryid) %>%
  mutate(lag_cases             = lag(total_cases, n=16),
         lag_cases             = ifelse(is.na(lag_cases)==TRUE, 0, lag_cases)) %>% 
  summarise(continent          = first(continent),
            location           = first(location),
            cases              = max(total_cases, na.rm=TRUE),
            lag_cases          = max(lag_cases, na.rm=TRUE),
            deaths             = max(total_deaths, na.rm=TRUE),
            tests              = max(total_tests, na.rm=TRUE),
            vaccinations       = max(total_vaccinations, na.rm=TRUE),
            stringency         = mean(s_stringency_index, na.rm=TRUE),
            reproduction       = mean(s_reproduction_rate, na.rm=TRUE),
            rate_at_lockdown   = first(rate_at_lockdown),
            time_to_lockdown   = first(time_to_lockdown),
            population         = first(population),
            gdp_per_capita     = first(gdp_per_capita),
            hospital_beds      = first(hospital_beds_per_thousand), 
            hdi                = first(human_development_index)) %>% 
  mutate(tests                 = ifelse(is.infinite(tests)==TRUE, NA, tests),
         vaccinations          = ifelse(is.infinite(vaccinations)==TRUE, NA, vaccinations))

COVID_data_Aug <- COVID_data_Aug %>% 
  group_by(countryid) %>%
  mutate(lag_cases             = lag(total_cases, n=16),
         lag_cases             = ifelse(is.na(lag_cases)==TRUE, 0, lag_cases)) %>% 
  summarise(continent          = first(continent),
            location           = first(location),
            cases              = max(total_cases, na.rm=TRUE),
            lag_cases          = max(lag_cases, na.rm=TRUE),
            deaths             = max(total_deaths, na.rm=TRUE),
            tests              = max(total_tests, na.rm=TRUE),
            vaccinations       = max(total_vaccinations, na.rm=TRUE),
            stringency         = mean(s_stringency_index, na.rm=TRUE),
            reproduction       = mean(s_reproduction_rate, na.rm=TRUE),
            rate_at_lockdown   = first(rate_at_lockdown),
            time_to_lockdown   = first(time_to_lockdown),
            population         = first(population),
            gdp_per_capita     = first(gdp_per_capita),
            hospital_beds      = first(hospital_beds_per_thousand), 
            hdi                = first(human_development_index)) %>% 
  mutate(tests                 = ifelse(is.infinite(tests)==TRUE, NA, tests),
         vaccinations          = ifelse(is.infinite(vaccinations)==TRUE, NA, vaccinations))

#Import 'exposure' data from the Global Burden of Disease 2019
NCD_data <- read_csv("C:/Users/Admin/Dropbox/Causal Insights/EIC/COVID Data/Data/GBDdownload.csv", show_col_types = FALSE)

#Drop superfluous collumns:
NCD_data        <- NCD_data %>% 
  select(c(measure_id, location_name, val))


#Convert into multiple columns per country 
NCD_data        <- NCD_data %>% 
  unite(ncds, measure_id) %>% 
  spread(ncds, val)

#Rename remaining variables
NCD_data        <- rename(NCD_data, location = location_name,
                          ncd_deaths = "1",
                          ncd_dalys = "2")

#Reduce to relevant countries
List_location  <- as.list(COVID_data$location) 

NCD_data  <- NCD_data %>%
  mutate(location = ifelse(location=="Bolivia (Plurinational State of)", "Bolivia", location),
         location = ifelse(location=="Côte d'Ivoire", "Cote d'Ivoire", location),
         location = ifelse(location=="Iran (Islamic Republic of)", "Iran", location),
         location = ifelse(location=="Russian Federation", "Russia", location),
         location = ifelse(location=="Syrian Arab Republic", "Syria", location),
         location = ifelse(location=="United Republic of Tanzania", "Tanzania", location),
         location = ifelse(location=="Venezuela (Bolivarian Republic of)", "Venezuela", location),
         location = ifelse(location=="Viet Nam", "Vietnam", location))
         
NCD_data        <- subset(NCD_data, location %in% List_location)

#Import other country-level data and reduce to relevant
Country_Info_data <- read_csv("C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Country_info.csv", show_col_types = FALSE)
Country_Info_data <- subset(Country_Info_data, countryid %in% List_countryid)

#Import age and sex data
Age_Sex_data <- read_csv("C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\age_and_sex.csv", show_col_types = FALSE)
Age_Sex_data <- subset(Age_Sex_data, countryid %in% List_countryid)

#Drop with missing age/sex data:
Age_Sex_data <- Age_Sex_data %>% 
  filter(is.na(m80)==FALSE)

#Join datasets
COVID_NCD_data  <- inner_join(COVID_data, NCD_data, by=c("location"))
COVID_NCD_data  <- inner_join(COVID_NCD_data, Country_Info_data, by=c("countryid"))
COVID_NCD_data  <- inner_join(COVID_NCD_data, Age_Sex_data, by=c("countryid"))

#And for sensitivity analyses
COVID_NCD_data_Jan  <- inner_join(COVID_data_Jan, NCD_data, by=c("location"))
COVID_NCD_data_Jan  <- inner_join(COVID_NCD_data_Jan, Country_Info_data, by=c("countryid"))
COVID_NCD_data_Jan  <- inner_join(COVID_NCD_data_Jan, Age_Sex_data, by=c("countryid"))

COVID_NCD_data_Aug  <- inner_join(COVID_data_Aug, NCD_data, by=c("location"))
COVID_NCD_data_Aug  <- inner_join(COVID_NCD_data_Aug, Country_Info_data, by=c("countryid"))
COVID_NCD_data_Aug  <- inner_join(COVID_NCD_data_Aug, Age_Sex_data, by=c("countryid"))

#Make new sex and age variables:
COVID_NCD_data    <- COVID_NCD_data %>% 
  mutate(m0_49    = (m0_14 + m15_19 + m20_24 + m25_29 + m30_34 + m35_39 + m40_44 + m45_49),
         m50_69   = (m50_54  + m55_59 + m60_64 + m65_69),
         m70      = (m70_74 + m75_79 + m80),
         f0_49    = (f0_14 + f15_19 + f20_24 + f25_29 + f30_34 + f35_39 + f40_44 + f45_49),
         f50_69   = (f50_54  + f55_59 + f60_64 + f65_69),
         f70      = (f70_74 + f75_79 + f80)) 

#Repeat in sensitivity analyses datasets:
COVID_NCD_data_Jan    <- COVID_NCD_data_Jan %>% 
  mutate(m0_49    = (m0_14 + m15_19 + m20_24 + m25_29 + m30_34 + m35_39 + m40_44 + m45_49),
         m50_69   = (m50_54  + m55_59 + m60_64 + m65_69),
         m70      = (m70_74 + m75_79 + m80),
         f0_49    = (f0_14 + f15_19 + f20_24 + f25_29 + f30_34 + f35_39 + f40_44 + f45_49),
         f50_69   = (f50_54  + f55_59 + f60_64 + f65_69),
         f70      = (f70_74 + f75_79 + f80)) 

COVID_NCD_data_Aug    <- COVID_NCD_data_Aug %>% 
  mutate(m0_49    = (m0_14 + m15_19 + m20_24 + m25_29 + m30_34 + m35_39 + m40_44 + m45_49),
         m50_69   = (m50_54  + m55_59 + m60_64 + m65_69),
         m70      = (m70_74 + m75_79 + m80),
         f0_49    = (f0_14 + f15_19 + f20_24 + f25_29 + f30_34 + f35_39 + f40_44 + f45_49),
         f50_69   = (f50_54  + f55_59 + f60_64 + f65_69),
         f70      = (f70_74 + f75_79 + f80)) 


#### DESCRIPTIVE ANALYSIS ####

#Save some important summary numbers:  

mean_population                <- mean(COVID_NCD_data$pop_total)
mean_cases                     <- mean(COVID_NCD_data$cases)
mean_lag_cases                 <- mean(COVID_NCD_data$lag_cases)
mean_deaths                    <- mean(COVID_NCD_data$deaths)
var_deaths                     <- var(COVID_NCD_data$deaths)
mean_ncd_deaths                <- mean(COVID_NCD_data$ncd_deaths)
mean_ncd_deaths_pop            <- mean_ncd_deaths/mean_population
mean_ncd_dalys                 <- mean(COVID_NCD_data$ncd_dalys)
mean_ncd_dalys_pop             <- mean_ncd_dalys/mean_population
mean_cfr                       <- mean(COVID_NCD_data$deaths)/mean(COVID_NCD_data$lag_cases)
variance_cfr                   <- var(COVID_NCD_data$deaths/COVID_NCD_data$lag_cases)

#Describe the sample

#Create variables for summarising
COVID_NCD_data <- COVID_NCD_data %>% 
  mutate(ncd_deaths_1k       = ncd_deaths/1000,
         ncd_mortality       = ncd_deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         ncd_dalys_1k        = ncd_dalys/1000,
         ncd_daly_ratio      = ncd_dalys/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_cases_1k      = cases/1000,
         covid_incidence     = cases/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_deaths_1k     = deaths/1000,
         covid_mortality     = deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*100000,
         covid_cfr           = (deaths/lag_cases)*100,
         total_pop_1k        = (m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)/1000,
         total_males         = (m0_49 + m50_69 + m70)/1000/total_pop_1k*100,
         total_0_49          = (m0_49 + f0_49)/1000/total_pop_1k*100,
         total_50_69         = (m50_69 + f50_69)/1000/total_pop_1k*100,
         total_70            = (m70 + f70)/1000/total_pop_1k*100,
         gdp_1k              = gdp_per_capita/1000,
         doctors_1k          = doctors_10k/10,
         tests_1m            = tests/1000000,
         country_area_1k     = country_area/1000,
         id                  = row_number())

#Repeat for sensitivity analysis
COVID_NCD_data_Jan <- COVID_NCD_data_Jan %>% 
  mutate(ncd_deaths_1k       = ncd_deaths/1000,
         ncd_mortality       = ncd_deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         ncd_dalys_1k        = ncd_dalys/1000,
         ncd_daly_ratio      = ncd_dalys/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_cases_1k      = cases/1000,
         covid_incidence     = cases/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_deaths_1k     = deaths/1000,
         covid_mortality     = deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*100000,
         covid_cfr           = (deaths/lag_cases)*100,
         total_pop_1k        = (m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)/1000,
         total_males         = (m0_49 + m50_69 + m70)/1000/total_pop_1k*100,
         total_0_49          = (m0_49 + f0_49)/1000/total_pop_1k*100,
         total_50_69         = (m50_69 + f50_69)/1000/total_pop_1k*100,
         total_70            = (m70 + f70)/1000/total_pop_1k*100,
         gdp_1k              = gdp_per_capita/1000,
         doctors_1k          = doctors_10k/10,
         tests_1m            = tests/1000000,
         country_area_1k     = country_area/1000,
         id                  = row_number())

COVID_NCD_data_Jan <- COVID_NCD_data_Jan[which(is.finite(COVID_NCD_data_Jan$deaths)==TRUE),] 

COVID_NCD_data_Aug <- COVID_NCD_data_Aug %>% 
  mutate(ncd_deaths_1k       = ncd_deaths/1000,
         ncd_mortality       = ncd_deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         ncd_dalys_1k        = ncd_dalys/1000,
         ncd_daly_ratio      = ncd_dalys/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_cases_1k      = cases/1000,
         covid_incidence     = cases/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_deaths_1k     = deaths/1000,
         covid_mortality     = deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*100000,
         covid_cfr           = (deaths/lag_cases)*100,
         total_pop_1k        = (m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)/1000,
         total_males         = (m0_49 + m50_69 + m70)/1000/total_pop_1k*100,
         total_0_49          = (m0_49 + f0_49)/1000/total_pop_1k*100,
         total_50_69         = (m50_69 + f50_69)/1000/total_pop_1k*100,
         total_70            = (m70 + f70)/1000/total_pop_1k*100,
         gdp_1k              = gdp_per_capita/1000,
         doctors_1k          = doctors_10k/10,
         tests_1m            = tests/1000000,
         country_area_1k     = country_area/1000,
         id                  = row_number())

# Restict to sample with complete data on core variables
COVID_NCD_data <- COVID_NCD_data[which(is.na(COVID_NCD_data$deaths)==FALSE &
                                       is.na(COVID_NCD_data$lag_cases)==FALSE &
                                       is.na(COVID_NCD_data$ncd_deaths)==FALSE & 
                                       is.na(COVID_NCD_data$ncd_dalys)==FALSE & 
                                       is.na(COVID_NCD_data$m0_49)==FALSE & 
                                       is.na(COVID_NCD_data$m50_69)==FALSE & 
                                       is.na(COVID_NCD_data$m70)==FALSE & 
                                       is.na(COVID_NCD_data$f0_49)==FALSE & 
                                       is.na(COVID_NCD_data$f50_69)==FALSE &
                                       is.na(COVID_NCD_data$f70)==FALSE),]

#List of variables to be summarised
List_Cont_Vars         <- list("ncd_deaths_1k", "ncd_mortality", "ncd_dalys_1k", "ncd_daly_ratio", "deaths", "covid_mortality", "covid_cases_1k", "covid_incidence", "covid_cfr", "total_pop_1k", "total_males", "total_0_49", "total_50_69", "total_70", "gdp_1k", "percent_urban", "smoking2018", "obesity2016", "doctors_1k", "hospital_beds", "healthcare_exp", "ph_exp", "stringency", "reproduction", "rate_at_lockdown", "time_to_lockdown", "tests_1m", "country_area_1k", "pollution")

#Create black dataframe to store summary information
Cont_Sum    <- data.frame(variable=character(0), min=numeric(0), max=numeric(0), median=numeric(0), q1=numeric(0), q3=numeric(0), n=numeric(0), stringsAsFactors=F)

for (i in 1:length(List_Cont_Vars)) {
  
  matrix        <-  cbind(NULL, List_Cont_Vars[i], round(as.numeric(length(COVID_NCD_data$id[which(is.na(COVID_NCD_data[,unlist(List_Cont_Vars[i])])==FALSE)])), digits=1), round(as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0, na.rm=TRUE)), digits=1), round(as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0.25, na.rm=TRUE)), digits=1), round(as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0.5, na.rm=TRUE)), digits=1), round(as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0.75, na.rm=TRUE)), digits=1), round(as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=1, na.rm=TRUE)), digits=1))
  Cont_Sum      <-  rbind(Cont_Sum, as.data.frame(matrix[1,], col.names = c("Variable", "N", "Min", "Q1", "Median", "Q3", "Max")))
}

write.csv(Cont_Sum, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Summary.csv")

#### PRIMARY ANALYSIS ####

#Now to determine the best adjustment set:
COVID_NCD_data$ncd_deaths <- round(COVID_NCD_data$ncd_deaths, digits=0)
COVID_NCD_data$ncd_dalys  <- round(COVID_NCD_data$ncd_dalys, digits=0)  

#Run algorthimic search for best model:
best_model <- glmulti(deaths ~  m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + gdp_per_capita + hospital_beds + smoking2018 + obesity2016 + doctors_10k + healthcare_exp + ph_exp + percent_urban + pollution, data = COVID_NCD_data, offset=log(COVID_NCD_data$cases), family=gaussian(link="log"), crit=bic, fitfunction=glm, method="l", level=1)
summary(best_model)

# The best model - according to the BIC - includes the age/sex variables plus smoking2018 and healthcare_exp




#Run estimation models



#MODEL 0 - No adjustment - crude associations


#Deaths model 0
model0_deaths                  <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k, data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model0_deaths)

#Call coefficient and CIs
exp(summary(model0_deaths)$coefficients[2,])
exp(confint(model0_deaths, "ncd_deaths_1k"))
#Transformed for interpretation:
exp(summary(model0_deaths)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model0_deaths, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)


#DALYs model 0
model0_dalys                  <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model0_dalys)

#Call coefficient and CIs
exp(summary(model0_dalys)$coefficients[2,])
exp(confint(model0_dalys, "I(ncd_dalys_1k/100)"))
#Transformed for interpretation:
exp(summary(model0_dalys)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model0_dalys, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)


#Repeat for sensitivity analysis
model0_deaths_Jan              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k, data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model0_deaths_Aug              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k, data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model0_deaths_Jan)$coefficients[2,])
exp(confint(model0_deaths_Jan, "ncd_deaths_1k"))
exp(summary(model0_deaths_Jan)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model0_deaths_Jan, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)
exp(summary(model0_deaths_Aug)$coefficients[2,])
exp(confint(model0_deaths_Aug, "ncd_deaths_1k"))
exp(summary(model0_deaths_Aug)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model0_deaths_Aug, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

model0_dalys_Jan               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model0_dalys_Aug               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model0_dalys_Jan)$coefficients[2,])
exp(confint(model0_dalys_Jan, "I(ncd_dalys_1k/100)"))
exp(summary(model0_dalys_Jan)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model0_dalys_Jan, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)
exp(summary(model0_dalys_Aug)$coefficients[2,])
exp(confint(model0_dalys_Aug, "I(ncd_dalys_1k/100)"))
exp(summary(model0_dalys_Aug)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model0_dalys_Aug, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)

#MODEL 1 - Age and sex adjustment only


#Deaths model 1
model1_deaths                  <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model1_deaths)
results1_deaths_rr             <- as.matrix(exp(summary(model1_deaths)$coefficients[,1]))
results1_deaths_lci            <- as.matrix(exp(confint(model1_deaths)[,1]))
results1_deaths_uci            <- as.matrix(exp(confint(model1_deaths)[,2]))

results1_deaths                <- cbind(results1_deaths_rr, results1_deaths_lci, results1_deaths_uci)

write.csv(results1_deaths, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model1Deaths.csv")

#Call coefficient and CIs
exp(summary(model1_deaths)$coefficients[2,])
exp(confint(model1_deaths, "ncd_deaths_1k"))
#Transformed for interpretation:
exp(summary(model1_deaths)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model1_deaths, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)


#DALYs model 1
model1_dalys                  <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model1_dalys)
results1_dalys_rr             <- as.matrix(exp(summary(model1_dalys)$coefficients[,1]))
results1_dalys_lci            <- as.matrix(exp(confint(model1_dalys)[,1]))
results1_dalys_uci            <- as.matrix(exp(confint(model1_dalys)[,2]))

results1_dalys                <- cbind(results1_dalys_rr, results1_dalys_lci, results1_dalys_uci)

write.csv(results1_dalys, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model1Dalys.csv")

#Call coefficient and CIs
exp(summary(model1_dalys)$coefficients[2,])
exp(confint(model1_dalys, "I(ncd_dalys_1k/100)"))
#Transformed for interpretation:
exp(summary(model1_dalys)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model1_dalys, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)



#Repeat for sensitivity analysis
model1_deaths_Jan              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model1_deaths_Aug              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model1_deaths_Jan)$coefficients[2,])
exp(confint(model1_deaths_Jan, "ncd_deaths_1k"))
exp(summary(model1_deaths_Jan)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model1_deaths_Jan, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)
exp(summary(model1_deaths_Aug)$coefficients[2,])
exp(confint(model1_deaths_Aug, "ncd_deaths_1k"))
exp(summary(model1_deaths_Aug)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model1_deaths_Aug, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

model1_dalys_Jan               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model1_dalys_Aug               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model1_dalys_Jan)$coefficients[2,])
exp(confint(model1_dalys_Jan, "I(ncd_dalys_1k/100)"))
exp(summary(model1_dalys_Jan)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model1_dalys_Jan, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)
exp(summary(model1_dalys_Aug)$coefficients[2,])
exp(confint(model1_dalys_Aug, "I(ncd_dalys_1k/100)"))
exp(summary(model1_dalys_Aug)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model1_dalys_Aug, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)


#MODEL 2 - Age, sex, and best confounders



#Deaths model 2
model2_deaths                  <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model2_deaths)
results2_deaths_rr             <- as.matrix(exp(summary(model2_deaths)$coefficients[,1]))
results2_deaths_lci            <- as.matrix(exp(confint(model2_deaths)[,1]))
results2_deaths_uci            <- as.matrix(exp(confint(model2_deaths)[,2]))

results2_deaths                <- cbind(results2_deaths_rr, results2_deaths_lci, results2_deaths_uci)

write.csv(results2_deaths, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model2Deaths.csv")

#Call coefficient and CIs
exp(summary(model2_deaths)$coefficients[2,])
exp(confint(model2_deaths, "ncd_deaths_1k"))
#Transformed for interpretation:
exp(summary(model2_deaths)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2_deaths, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)


#DALYs model 2
model2_dalys                  <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model2_dalys)
results2_dalys_rr             <- as.matrix(exp(summary(model2_dalys)$coefficients[,1]))
results2_dalys_lci            <- as.matrix(exp(confint(model2_dalys)[,1]))
results2_dalys_uci            <- as.matrix(exp(confint(model2_dalys)[,2]))

results2_dalys                <- cbind(results2_dalys_rr, results2_dalys_lci, results2_dalys_uci)

write.csv(results2_dalys, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model2Dalys.csv")

#Call coefficient and CIs
exp(summary(model2_dalys)$coefficients[2,])
exp(confint(model2_dalys, "I(ncd_dalys_1k/100)"))
#Transformed for interpretation:
exp(summary(model2_dalys)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2_dalys, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)



#Repeat for sensitivity analysis
model2_deaths_Jan              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model2_deaths_Aug              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model2_deaths_Jan)$coefficients[2,])
exp(confint(model2_deaths_Jan, "ncd_deaths_1k"))
exp(summary(model2_deaths_Jan)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2_deaths_Jan, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)
exp(summary(model2_deaths_Aug)$coefficients[2,])
exp(confint(model2_deaths_Aug, "ncd_deaths_1k"))
exp(summary(model2_deaths_Aug)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2_deaths_Aug, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

model2_dalys_Jan               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model2_dalys_Aug               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model2_dalys_Jan)$coefficients[2,])
exp(confint(model2_dalys_Jan, "I(ncd_dalys_1k/100)"))
exp(summary(model2_dalys_Jan)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2_dalys_Jan, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)
exp(summary(model2_dalys_Aug)$coefficients[2,])
exp(confint(model2_dalys_Aug, "I(ncd_dalys_1k/100)"))
exp(summary(model2_dalys_Aug)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2_dalys_Aug, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)




#MODELS 2A+B - Age, sex, best confounders, and testing

#Deaths model 2A - sample with testing data
model2a_deaths              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data[which(is.na(COVID_NCD_data$tests)==FALSE),], family=quasipoisson(), maxit = 100)
summary(model2a_deaths)
results2a_deaths_rr             <- as.matrix(exp(summary(model2a_deaths)$coefficients[,1]))
results2a_deaths_lci            <- as.matrix(exp(confint(model2a_deaths)[,1]))
results2a_deaths_uci            <- as.matrix(exp(confint(model2a_deaths)[,2]))

results2a_deaths                <- cbind(results2a_deaths_rr, results2a_deaths_lci, results2a_deaths_uci)

write.csv(results2a_deaths, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model2ADeaths.csv")

#Call coefficient and CIs
exp(model2a_deaths$coefficients[2])
exp(confint(model2a_deaths, "ncd_deaths_1k"))
#Transformed for interpretation:
exp(summary(model2a_deaths)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2a_deaths, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)



#DALYs model 2A - sample with testing data
model2a_dalys                  <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data[which(is.na(COVID_NCD_data$tests)==FALSE),], family=quasipoisson(), maxit = 100)
summary(model2a_dalys)
results2a_dalys_rr             <- as.matrix(exp(summary(model2a_dalys)$coefficients[,1]))
results2a_dalys_lci            <- as.matrix(exp(confint(model2a_dalys)[,1]))
results2a_dalys_uci            <- as.matrix(exp(confint(model2a_dalys)[,2]))

results2a_dalys                <- cbind(results2a_dalys_rr, results2a_dalys_lci, results2a_dalys_uci)

write.csv(results2a_dalys, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model2ADalys.csv")

#Call coefficient and CIs
exp(model2a_dalys$coefficients[2])
exp(confint(model2a_dalys, "I(ncd_dalys_1k/100)"))
#Transformed for interpretation:
exp(summary(model2a_dalys)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2a_dalys, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)


#Repeat for sensitivity analysis
model2a_deaths_Jan              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Jan[which(is.na(COVID_NCD_data_Jan$tests)==FALSE),], family=quasipoisson(), maxit = 100)
model2a_deaths_Aug              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Aug[which(is.na(COVID_NCD_data_Aug$tests)==FALSE),], family=quasipoisson(), maxit = 100)

exp(summary(model2a_deaths_Jan)$coefficients[2,])
exp(confint(model2a_deaths_Jan, "ncd_deaths_1k"))
exp(summary(model2a_deaths_Jan)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2a_deaths_Jan, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)
exp(summary(model2a_deaths_Aug)$coefficients[2,])
exp(confint(model2a_deaths_Aug, "ncd_deaths_1k"))
exp(summary(model2a_deaths_Aug)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2a_deaths_Aug, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

model2a_dalys_Jan               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Jan[which(is.na(COVID_NCD_data_Jan$tests)==FALSE),], family=quasipoisson(), maxit = 100)
model2a_dalys_Aug               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000), data=COVID_NCD_data_Aug[which(is.na(COVID_NCD_data_Aug$tests)==FALSE),], family=quasipoisson(), maxit = 100)

exp(summary(model2a_dalys_Jan)$coefficients[2,])
exp(confint(model2a_dalys_Jan, "I(ncd_dalys_1k/100)"))
exp(summary(model2a_dalys_Jan)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2a_dalys_Jan, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)
exp(summary(model2a_dalys_Aug)$coefficients[2,])
exp(confint(model2a_dalys_Aug, "I(ncd_dalys_1k/100)"))
exp(summary(model2a_dalys_Aug)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2a_dalys_Aug, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)



#Deaths model 2B - including testing data
model2b_deaths              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000) + tests_1m, data=COVID_NCD_data[which(is.na(COVID_NCD_data$tests)==FALSE),], family=quasipoisson(), maxit = 100)
summary(model2b_deaths)
results2b_deaths_rr             <- as.matrix(exp(summary(model2b_deaths)$coefficients[,1]))
results2b_deaths_lci            <- as.matrix(exp(confint(model2b_deaths)[,1]))
results2b_deaths_uci            <- as.matrix(exp(confint(model2b_deaths)[,2]))

results2b_deaths                <- cbind(results2b_deaths_rr, results2b_deaths_lci, results2b_deaths_uci)

write.csv(results2b_deaths, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model2BDeaths.csv")

#Call coefficient and CIs
exp(model2b_deaths$coefficients[2])
exp(confint(model2b_deaths, "ncd_deaths_1k"))
#Transformed for interpretation:
exp(summary(model2b_deaths)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2b_deaths, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)


#DALYs model 2B - including testing data
model2b_dalys              <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000) + tests_1m, data=COVID_NCD_data[which(is.na(COVID_NCD_data$tests)==FALSE),], family=quasipoisson(), maxit = 100)
summary(model2b_dalys)
results2b_dalys_rr             <- as.matrix(exp(summary(model2b_dalys)$coefficients[,1]))
results2b_dalys_lci            <- as.matrix(exp(confint(model2b_dalys)[,1]))
results2b_dalys_uci            <- as.matrix(exp(confint(model2b_dalys)[,2]))

results2b_dalys                <- cbind(results2b_dalys_rr, results2b_dalys_lci, results2b_dalys_uci)

write.csv(results2b_dalys, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\Model2BDalys.csv")

#Call coefficient and CIs
exp(model2b_dalys$coefficients[2])
exp(confint(model2b_dalys, "I(ncd_dalys_1k/100)"))
#Transformed for interpretation:
exp(summary(model2b_dalys)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2b_dalys, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)


#Repeat for sensitivity analysis
model2b_deaths_Jan              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000) + tests_1m, data=COVID_NCD_data_Jan[which(is.na(COVID_NCD_data_Jan$tests)==FALSE),], family=quasipoisson(), maxit = 100)
model2b_deaths_Aug              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000) + tests_1m, data=COVID_NCD_data_Aug[which(is.na(COVID_NCD_data_Aug$tests)==FALSE),], family=quasipoisson(), maxit = 100)

exp(summary(model2b_deaths_Jan)$coefficients[2,])
exp(confint(model2b_deaths_Jan, "ncd_deaths_1k"))
exp(summary(model2b_deaths_Jan)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2b_deaths_Jan, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)
exp(summary(model2b_deaths_Aug)$coefficients[2,])
exp(confint(model2b_deaths_Aug, "ncd_deaths_1k"))
exp(summary(model2b_deaths_Aug)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model2b_deaths_Aug, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

model2b_dalys_Jan               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000) + tests_1m, data=COVID_NCD_data_Jan[which(is.na(COVID_NCD_data_Jan$tests)==FALSE),], family=quasipoisson(), maxit = 100)
model2b_dalys_Aug               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + smoking2018 + I(healthcare_exp/1000) + tests_1m, data=COVID_NCD_data_Aug[which(is.na(COVID_NCD_data_Aug$tests)==FALSE),], family=quasipoisson(), maxit = 100)

exp(summary(model2b_dalys_Jan)$coefficients[2,])
exp(confint(model2b_dalys_Jan, "I(ncd_dalys_1k/100)"))
exp(summary(model2b_dalys_Jan)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2b_dalys_Jan, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)
exp(summary(model2b_dalys_Aug)$coefficients[2,])
exp(confint(model2b_dalys_Aug, "I(ncd_dalys_1k/100)"))
exp(summary(model2b_dalys_Aug)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model2b_dalys_Aug, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)

#MODEL 3 - Adjusting for all confounders and competing exposures 

#Deaths model 3

model3_deaths              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + gdp_1k + percent_urban + smoking2018 + obesity2016 + doctors_1k + hospital_beds + I(healthcare_exp/1000) + ph_exp + I(reproduction*10) + stringency + I(rate_at_lockdown/10) + I(time_to_lockdown/28) + pollution + I(country_area_1k/1000), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model3_deaths)
results3_deaths_rr             <- as.matrix(exp(summary(model3_deaths)$coefficients[,1]))
results3_deaths_lci            <- as.matrix(exp(confint(model3_deaths)[,1]))
results3_deaths_uci            <- as.matrix(exp(confint(model3_deaths)[,2]))

results3_deaths                <- cbind(results3_deaths_rr, results3_deaths_lci, results3_deaths_uci)

write.csv(results3_deaths, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\model3Deaths.csv")

#Call coefficient and CIs
exp(model3_deaths$coefficients[2])
exp(confint(model3_deaths, "ncd_deaths_1k"))

#Transformed for interpretation:
exp(summary(model3_deaths)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model3_deaths, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

#DALYs model 3
model3_dalys              <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + gdp_1k + percent_urban + smoking2018 + obesity2016 + doctors_1k + hospital_beds + I(healthcare_exp/1000) + ph_exp + I(reproduction*10) + stringency + I(rate_at_lockdown/10) + I(time_to_lockdown/28) + pollution + I(country_area_1k/1000), data=COVID_NCD_data, family=quasipoisson(), maxit = 100)
summary(model3_dalys)
results3_dalys_rr             <- as.matrix(exp(summary(model3_dalys)$coefficients[,1]))
results3_dalys_lci            <- as.matrix(exp(confint(model3_dalys)[,1]))
results3_dalys_uci            <- as.matrix(exp(confint(model3_dalys)[,2]))

results3_dalys                <- cbind(results3_dalys_rr, results3_dalys_lci, results3_dalys_uci)
write.csv(results3_dalys, file = "C:\\Users\\Admin\\Dropbox\\Causal Insights\\EIC\\COVID Data\\Data\\model3Dalys.csv")

#Call coefficient and CIs
exp(model3_dalys$coefficients[2])
exp(confint(model3_dalys, "I(ncd_dalys_1k/100)"))

#Transformed for interpretation:
exp(summary(model3_dalys)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model3_dalys, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)


#Repeat for sensitivity analysis
model3_deaths_Jan              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + gdp_1k + percent_urban + smoking2018 + obesity2016 + doctors_1k + hospital_beds + I(healthcare_exp/1000) + ph_exp + I(reproduction*10) + stringency + I(rate_at_lockdown/10) + I(time_to_lockdown/28) + pollution + I(country_area_1k/1000), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model3_deaths_Aug              <- glm(deaths ~ offset(log(lag_cases)) + ncd_deaths_1k + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + gdp_1k + percent_urban + smoking2018 + obesity2016 + doctors_1k + hospital_beds + I(healthcare_exp/1000) + ph_exp + I(reproduction*10) + stringency + I(rate_at_lockdown/10) + I(time_to_lockdown/28) + pollution + I(country_area_1k/1000), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model3_deaths_Jan)$coefficients[2,])
exp(confint(model3_deaths_Jan, "ncd_deaths_1k"))
exp(summary(model3_deaths_Jan)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model3_deaths_Jan, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)
exp(summary(model3_deaths_Aug)$coefficients[2,])
exp(confint(model3_deaths_Aug, "ncd_deaths_1k"))
exp(summary(model3_deaths_Aug)$coefficients[2]*0.1*mean_ncd_deaths/1000)
exp(confint(model3_deaths_Aug, "ncd_deaths_1k")*0.1*mean_ncd_deaths/1000)

model3_dalys_Jan               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + gdp_1k + percent_urban + smoking2018 + obesity2016 + doctors_1k + hospital_beds + I(healthcare_exp/1000) + ph_exp + I(reproduction*10) + stringency + I(rate_at_lockdown/10) + I(time_to_lockdown/28) + pollution + I(country_area_1k/1000), data=COVID_NCD_data_Jan, family=quasipoisson(), maxit = 100)
model3_dalys_Aug               <- glm(deaths ~ offset(log(lag_cases)) + I(ncd_dalys_1k/100) + I(m0_49/100000) + I(m50_69/100000) + I(m70/100000) + I(f0_49/100000) + I(f50_69/100000) + I(f70/100000) + gdp_1k + percent_urban + smoking2018 + obesity2016 + doctors_1k + hospital_beds + I(healthcare_exp/1000) + ph_exp + I(reproduction*10) + stringency + I(rate_at_lockdown/10) + I(time_to_lockdown/28) + pollution + I(country_area_1k/1000), data=COVID_NCD_data_Aug, family=quasipoisson(), maxit = 100)

exp(summary(model3_dalys_Jan)$coefficients[2,])
exp(confint(model3_dalys_Jan, "I(ncd_dalys_1k/100)"))
exp(summary(model3_dalys_Jan)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model3_dalys_Jan, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)
exp(summary(model3_dalys_Aug)$coefficients[2,])
exp(confint(model3_dalys_Aug, "I(ncd_dalys_1k/100)"))
exp(summary(model3_dalys_Aug)$coefficients[2]*0.1*mean_ncd_dalys/100000)
exp(confint(model3_dalys_Aug, "I(ncd_dalys_1k/100)")*0.1*mean_ncd_dalys/100000)


#### GRAPH ####

#Determine deaths/CFR predicted by age, sex, smoking, and healthcare_exp
model_predicted_deaths             <- glm(deaths ~ offset(log(lag_cases)) + I(m0_49-mean(m0_49)) + I(m50_69-mean(m50_69)) + I(m70-mean(m70)) + I(f0_49-mean(f0_49)) + I(f50_69-mean(f50_69)) + I(f70-mean(f70)) + I(smoking2018-mean(smoking2018, na.rm=TRUE)) + I(healthcare_exp-mean(healthcare_exp, na.rm=TRUE)), data=COVID_NCD_data, family=gaussian(link=log), na.action=na.exclude, maxit=100)
COVID_NCD_data$pred_deaths         <- predict(model_predicted_deaths, type="response")
COVID_NCD_data$res_deaths          <- residuals(model_predicted_deaths, type="response")
COVID_NCD_data$res_cfr             <- (COVID_NCD_data$res_deaths/COVID_NCD_data$lag_cases)*100

#Now plot the relationship between NCD mortality/DALY ratio and cfr/res_cfr:

options(scipen=999)

#NCD Mortaltiy - CRUDE CFR

ggplot(data = COVID_NCD_data) +
  theme_classic(base_size = 12,
                base_family="sans") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(colour="black", size=14),
        axis.text.y=element_text(colour="black", size=14),
        axis.title.x =element_text(size=14, face="bold"),
        axis.title.y =element_text(size=14, face="bold")) +
  geom_smooth(aes(x = ncd_mortality, y = covid_cfr), colour="#C01E35", method="loess", span=1, fill="#FAE5E3", size=1) +
  geom_point(aes(x = ncd_mortality, y = covid_cfr), colour="#C01E35", size=0.5) +
  scale_x_continuous(name="NCD Mortality Ratio (per Thousand)", expand=c(0,0), n.breaks=5, limits=c(0, 20)) +
  scale_y_continuous(name="Observed COVID-19 CFR %", expand=c(0,0), n.breaks=5, limits=c(0, 10),  labels = scales::number_format(accuracy = 0.1))

#NCD DALY Ratio - CRUDE CFR

ggplot(data = COVID_NCD_data) +
  theme_classic(base_size = 12,
                base_family="sans") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(colour="black", size=14),
        axis.text.y=element_text(colour="black", size=14),
        axis.title.x =element_text(size=14, face="bold"),
        axis.title.y =element_text(size=14, face="bold")) +
  geom_smooth(aes(x = ncd_daly_ratio, y = covid_cfr), colour="#C01E35", method="loess", span=1, fill="#FAE5E3", size=1) +
  geom_point(aes(x = ncd_daly_ratio, y = covid_cfr), colour="#C01E35", size=0.5) +
  scale_x_continuous(name="NCD DALY Ratio (per Thousand)", expand=c(0,0), n.breaks=5, limits=c(100, 500)) +
  scale_y_continuous(name="Observed COVID-19 CFR %", expand=c(0,0), n.breaks=5, limits=c(0, 10), labels = scales::number_format(accuracy = 0.1))

#NCD Mortaltiy - Residual CFR

ggplot(data = COVID_NCD_data) +
  theme_classic(base_size = 12,
                base_family="sans") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(colour="black", size=14),
        axis.text.y=element_text(colour="black", size=14),
        axis.title.x =element_text(size=14, face="bold"),
        axis.title.y =element_text(size=14, face="bold")) +
  geom_smooth(aes(x = ncd_mortality, y = res_cfr), colour="#C01E35", method="loess", span=1, fill="#FAE5E3", size=1) +
  geom_point(aes(x = ncd_mortality, y = res_cfr), colour="#C01E35", size=0.5) +
  scale_x_continuous(name="NCD Mortality Ratio (per Thousand)", expand=c(0,0), n.breaks=5, limits=c(0, 20)) +
  scale_y_continuous(name="Residual COVID-19 CFR % (Observed-Expected)", expand=c(0,0), n.breaks=5, limits=c(-10, 10), labels = scales::number_format(accuracy = 0.1))

#NCD DALY ratio - Residual CFR

ggplot(data = COVID_NCD_data) +
  theme_classic(base_size = 12,
                base_family="sans") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(colour="black", size=14),
        axis.text.y=element_text(colour="black", size=14),
        axis.title.x =element_text(size=14, face="bold"),
        axis.title.y =element_text(size=14, face="bold")) +
  geom_smooth(aes(x = ncd_daly_ratio, y = res_cfr), colour="#C01E35", method="loess", span=1, fill="#FAE5E3", size=1) +
  geom_point(aes(x = ncd_daly_ratio, y = res_cfr), colour="#C01E35", size=0.5) +
  scale_x_continuous(name="NCD DALY Ratio (per Thousand)", expand=c(0,0), n.breaks=5, limits=c(100, 500)) +
  scale_y_continuous(name="Residual COVID-19 CFR % (Observed-Expected)", expand=c(0,0), n.breaks=5, limits=c(-10, 10), labels = scales::number_format(accuracy = 0.1))
