# Data set creation for NCWELL data
# For privacy reasons, please see links in manuscript for NCWELL database
# Murray et al data not public

rm(list=ls())
setwd("~/GitHub/Pairwise_Additive_GP_Tests")
# setwd("H:/Pairwise_Additive_GP_Tests")
library(tidyverse)
library(DescTools)

################################################################################
##################################### Demographic data #########################
################################################################################

tract_covar1 <- read.csv("NCWell/tract/tract_inputs_2019_v2.csv",
                         colClasses=list(Tract.2="character",
                                         County.FIPS="character"))

head(tract_covar1)
str(tract_covar1)

tract_covar2 <- tract_covar1 %>%
  rename(Age=Median.age..years., Black=Black.or.African.American, 
         Hispanic=Hispanic.or.Latino..of.any.race.,
         NativeAm=American.Indian.and.Alaska.Native) %>%
  mutate(across(all_of(c("Male", "Female", "Age",
                         "Black", "Hispanic", "Median.Income", "NativeAm")),
                as.numeric),
         Pop=Male+Female,
         Female.Pct=Female/Pop,
         Black.Pct=Black/Pop,
         Hisp.Pct=Hispanic/Pop,
         NatAm.Pct = NativeAm/Pop) 

head(tract_covar2)
str(tract_covar2)

convert_tract_fips <- function(tract.code){
  period <- grepl(".", tract.code, fixed=TRUE)
  if(!period){
    tract.code2 <- paste(tract.code, "00", sep="")
  }else{
    tract.code2 <- gsub(".", "", tract.code,fixed=TRUE)
  }

  tract.code3 <- tract.code2
  if(nchar(tract.code3)>6){
    stop("Long tract code")
  }
  while(nchar(tract.code3) != 6){
    tract.code3 <- paste("0", tract.code3, sep="")
  }
  return(tract.code3)
}


tract_covar2$tract.only.fips <- sapply(tract_covar2$Tract.2, FUN=convert_tract_fips)
tract_covar2$tract_fips <- paste(tract_covar2$County.FIPS, tract_covar2$tract.only.fips, sep="")


tract_covar <- tract_covar2 %>%
  select(tract_fips, Age, Median.Income, Female.Pct, 
         Black.Pct, Hisp.Pct, NatAm.Pct, Pop)



################################################################################
#                           Well use data                                      #
################################################################################

# Not public
well_use <- read.csv("NCWell//well_usage//Murray_output.csv") %>%
  mutate(tract_fips = as.character(geoid_10)) %>%
  select(tract_fips, pct_well=avg_pct_on_well_water)

################################################################################
#################### Health outcome (diabetes) and obesity  ####################
################################################################################

places_data <- read.csv("NCWell/PLACES__Local_Data_for_Better_Health,_Census_Tract_Data_2021_release_20260313.csv")

head(places_data)
diabetes_data <- places_data %>%
  filter(Measure=="Diagnosed diabetes among adults aged >=18 years") %>%
  select(tract_fips=LocationName, diabetes_pct=Data_Value, Geolocation) %>%
  mutate(tract_fips=as.character(tract_fips))
head(diabetes_data)

obesity_data <- places_data %>%
  filter(Measure=="Obesity among adults aged >=18 years") %>%
  select(tract_fips=LocationName, obesity_pct=Data_Value, Geolocation) %>%
  mutate(tract_fips=as.character(tract_fips))
head(obesity_data)


# diabetes_data %>%
#   distinct(Data_Value_Type)

str(diabetes_data)
extract_lat_lon <- function(geolocation){
  space.pos1 <- str_locate_all(pattern =' ',geolocation)
  space.pos <- space.pos1[[1]][2,1]
  
  left.par.pos1 <- str_locate_all(pattern ='\\(',geolocation)
  left.par.pos <- left.par.pos1[[1]][1,1]
  
  right.par.pos1 <- str_locate_all(pattern ='\\)',geolocation)
  right.par.pos <- right.par.pos1[[1]][1,1]
  
  # longitude is first number in these data
  lon <- substr(geolocation, start=left.par.pos + 1, space.pos - 1)
  lat <- substr(geolocation, start=space.pos+1, right.par.pos - 1)
  
  out <- c(lat, lon)
  
  out <- as.numeric(out)
  names(out) <- c("lat", "lon")
  return(out)
}
diabetes_data$lat <- rep(NA, nrow(diabetes_data))
diabetes_data$lon <- rep(NA, nrow(diabetes_data))

for(i in 1:nrow(diabetes_data)){
  latlon <- extract_lat_lon(diabetes_data$Geolocation[i])
  diabetes_data$lat[i] <- latlon["lat"]
  diabetes_data$lon[i] <- latlon["lon"]
}

diabetes_data %>%
  left_join(obesity_data) %>%
  drop_na() %>%
  nrow()

diabetes_data <- left_join(diabetes_data, obesity_data)


head(diabetes_data)

################################################################################
################ Exposure data (metals) ########################################
################################################################################

metals_dat_1 <- read.csv("NCWell/doi-10.15139-s3-bdqg9o/NCwelldata_distributions_tract_level.csv")


metals_dat_2 <- metals_dat_1 %>%
  select(tract_fips=FIPS, ends_with("Mean_avg")) %>%
  mutate(tract_fips=as.character(tract_fips))

head(metals_dat_2)
str(metals_dat_2)

nonzero.means <- colMeans(!is.na(metals_dat_2)) # proportion of tracts with wells?
detection.min <- (apply(metals_dat_2[,-1], FUN=min, na.rm=TRUE, MARGIN=2))
detection.rates <- rep(NA, length(detection.min))
names(detection.rates) <- names(detection.min)
metals_dat_3 <- metals_dat_2

metals_dat_3[is.na(metals_dat_3)] <- 0

for(j in 1:length(detection.min)){
  detection.rates[j] <- mean(metals_dat_3[,(j+1)] > detection.min[j])
}
detection.rates
sum(detection.rates > .3)
sort(detection.rates, decreasing = TRUE)

metals_dat <- metals_dat_2[,c(1,1+which(detection.rates > .7))]
head(metals_dat)
rowsums.metals <- rowSums(metals_dat[,2:ncol(metals_dat)])
metals_dat[rowsums.metals > 0,]

metals_dat %>%
  drop_na() %>%
  nrow()

## Create dataset filtered to wells with >= 4 tests (_res, restricted)
metals_dat_numtests_1 <- metals_dat_1[,c(1,1+which(detection.rates > .7))]  %>%
  mutate(tract_fips=as.character(FIPS))
head(metals_dat_numtests_1)
metals_dat_numtests_2 <- apply(sapply(metals_dat_numtests_1[,-1], as.numeric), MARGIN=1, FUN="min")
metals_dat_numtests <- data.frame(tract_fips=as.character(metals_dat_numtests_1[,1]), numtests=metals_dat_numtests_2)
head(metals_dat_numtests)
sum(metals_dat_numtests_2 > 3)

# add num tests to full metals data
metals_dat <- left_join(metals_dat, metals_dat_numtests)


hist(metals_dat$Arsenic.Mean_avg)
################################################################################
########################### Final data set #####################################
################################################################################

data_tract <- inner_join(diabetes_data, tract_covar) %>%
  inner_join(metals_dat) %>%
  inner_join(well_use) %>%
  select(-Geolocation) %>%
  drop_na()



write.csv(data_tract, file = "data/NCWell_tract.csv", row.names = FALSE)


