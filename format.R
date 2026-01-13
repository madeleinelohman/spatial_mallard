
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("beepr", "colorspace", "raster", "rstan","sf", "tidyverse")
ipak(packages)

### Set wd
setwd("~/Desktop/format_080126")
load('map_dat.RData')



#~~~~~~~~~~~~~~~~
# Define parameters of data
#~~~~~~~~~~~~~~~~
### Species
species <- "MALL"

### Years needed
start <- 1974
end <- 2024
official.end <- 2020
years <- start:end
n.years <- length(start:end)

#~~~~~~~~~~~~~~~~~~~~~~
# Spatial data
#~~~~~~~~~~~~~~~~~~~~~~
### Load in
shape <- grd
#st_crs(grd)
### Number of strata
n.strata <- nrow(shape)
### Determine projection
proj <- st_crs(shape)
shape$id <- 1:nrow(shape)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data cleaning and subsetting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- read.csv('NABBP_2025_grp_08.csv')
colnames(d) <- tolower(colnames(d))
d <- d %>%
  filter(!is.na(lat_dd) | !is.na(lon_dd))

#~~~~~~~~~~~~~~~~
# Appropriately subset data
#~~~~~~~~~~~~~~~~
### Create banding data
b <- d[which(d$record_source == "B"),]
### Create encounter data
# Obtained by shooting and reported by hunter
e <- d %>%
  filter(how_obtained_code == 1 & who_obtained_code == 21 & record_source == "E")


#~~~~~~~~~~~~~~~~~~~~~~
# Releases subsetting
#~~~~~~~~~~~~~~~~~~~~~~
### Banded during appropriate years
releases <- b %>%
  filter(event_year >= start & event_year < end)
  

### Released as a 'Normal wild bird'
releases <- releases %>%
  filter(bird_status == 3)

### Type of study captured for/how
# 00: Federal numbered metal band only.
# 04: Control band - For use in conjunction with reward band studies only. 
# 70: Spotlighted. 
capt.method <- c(0, 4, 70)
releases$extra_info_code <- as.numeric(releases$extra_info_code)
releases <- releases %>%
  filter(extra_info_code %in% capt.method)

### Type of band
#table(releases$BAND_TYPE)
# 01: aluminum\butt-end toll free
# 04: aluminum\butt-end new address
# 08: aluminum\pre-open toll free
# 11: aluminum\butt end
# 18: aluminum\pre-open
# 21: monel\butt end
# 23: monel\butt-end toll free
# 41: aluminum\butt-end (toll-free /web address)
# 51: incoloy or Stainless\butt-end
# 53: incoloy or Stainless\butt-end toll free
# W1: Aluminum butt-end web address

band.type <- c('01', '04', '08', '11', '18', '21', '23', '41', '51',
               '53', 'W1')
releases <- releases %>%
  filter(band_type_code %in% band.type)

### Only individuals banded during pre-season banding
# table(releases$EVENT_MONTH)
rel.months <- 7:9
releases <- releases %>%
  filter(event_month %in% rel.months)

#~~~~~~~~
# Spatially subset releases
#~~~~~~~~
### Convert to sf object
releases <- st_as_sf(releases, coords=c("lon_dd", "lat_dd"), crs=st_crs(4269),
                     remove=F)
### Convert to Albers projection
releases <- st_transform(releases, crs=proj)
# releases %>%
#   group_by(geometry) %>%
#   summarize(geo = unique(geometry)) %>%
#   ggplot() +
#   geom_sf()
### Find points that intersect with strata
rels <- st_join(shape, releases)
#rels <- st_join(shape, releases, join=st_contains)
# summary(rels$id)
rels <- rels[which(!is.na(rels$id)),]

#~~~~~~~~~~~~~~~~~~~~~~
# Encounters subsetting
#~~~~~~~~~~~~~~~~~~~~~~
### Get only bands from the releases dataset
encounters <- e %>%
  filter(band %in% rels$band)

### Only dead birds
cond <- c(3, 4, 5)
encounters <- encounters %>%
  filter(present_condition %in% cond)

### Only individuals shot during the actual months or actual hunting season
encounters <- encounters %>%
  filter(event_month <= 12 & (event_month < 4 | event_month > 8))

### Only individuals harvested before the end of the study
encounters <- encounters %>%
  filter(event_year <= end)

encounters$id <- st_drop_geometry(rels$id[match(encounters$band, rels$band)])
encounters$sex <- st_drop_geometry(rels$sex_code[match(encounters$band, rels$band)])

#~~~~~~~~
# Make encounters spatial
#~~~~~~~~
### As above
encounters <- st_as_sf(encounters, coords=c("lon_dd", "lat_dd"), crs=st_crs(4269))
encs <- st_transform(encounters, crs=proj)
# ggplot(encs) +
#   geom_sf()


### Get rid of unneeded objects
rm(list=c("capt.method", "band.type", "rel.months", "cond", 'b', 'e'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Form m-arrays
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~
# Subset into adults and juveniles
#~~~~~~~~~~~~~~~~~~~~~~
### Adults
a.ages <- c(1, 5, 6, 8)
rels.a <- rels %>%
  filter(age_code %in% a.ages)
encs.a <- encs %>%
  filter(band %in% rels.a$band)


### Juveniles
j.ages <- c(2, 4)
rels.j <- rels %>%
  filter(age_code %in% j.ages)
encs.j <- encs %>%
  filter(band %in% rels.j$band)


#~~~~~~~~~~~~~~~~~~~~~~
# M-arrays
#~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~
# Make function
#~~~~~~~~

marr <- function(rels, encs, n.years, n.strata, years){
  
  ### Sort birds into the year associated with the START of the hunting season they were shot in
  encs$event_year2 <- encs$event_year
  encs$event_year2[which(encs$event_month < 4)] <- encs$event_year[which(encs$event_month < 4)] - 1
  
  
  marr <- array(0, dim = c(n.years, n.years+1, n.strata)) # empty matrix
  
  for(s in 1:n.strata){
    for(i in 1:n.years){
      rels.want <- rels %>%
        filter(event_year == years[i] & id == s)
      
      for(t in i:n.years){
        encs.want <- encs %>%
          filter(band %in% rels.want$band & encs$event_year2 == years[t])
        
        marr[i,t,s] <- nrow(encs.want)
      }
      
      marr[i,n.years+1,s] <- nrow(rels.want) - sum(marr[i,,s])
    }
  }
  
  
  return(marr)
  
}


#~~~~~~~~
# Create m-arrays
#~~~~~~~~
marr.af <- marr(subset(rels.a, sex_code == 5), subset(encs.a, sex == 5), n.years, n.strata, years)
marr.jf <- marr(subset(rels.j, sex_code == 5), subset(encs.j, sex == 5), n.years, n.strata, years)
marr.am <- marr(subset(rels.a, sex_code == 4), subset(encs.a, sex == 4), n.years, n.strata, years)
marr.jm <- marr(subset(rels.j, sex_code == 4), subset(encs.j, sex == 4), n.years, n.strata, years)

#~~~~~~~~
# Sample sizes
#~~~~~~~~
### Females
## Total releases
sum(rowSums(marr.af))
sum(rowSums(marr.jf))
## Total encounters
sum(rowSums(marr.af[,1:(ncol(marr.af)-1),]))
sum(rowSums(marr.jf[,1:(ncol(marr.jf)-1),]))
# Proportion encounters
sum(rowSums(marr.af[,1:(ncol(marr.af)-1),])) / sum(rowSums(marr.af))
sum(rowSums(marr.jf[,1:(ncol(marr.jf)-1),])) / sum(rowSums(marr.jf))

### Males
## Total releases
sum(rowSums(marr.am))
sum(rowSums(marr.jm))
## Total encounters
sum(rowSums(marr.am[,1:(ncol(marr.am)-1),]))
sum(rowSums(marr.jm[,1:(ncol(marr.jm)-1),]))
# Proportion encounters
sum(rowSums(marr.af[,1:(ncol(marr.am)-1),])) / sum(rowSums(marr.am))
sum(rowSums(marr.jf[,1:(ncol(marr.jm)-1),])) / sum(rowSums(marr.jm))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(d, encounters, encs, encs.a, encs.j, proj, releases, rels, 
   rels.a, rels.j, a.ages, j.ages, packages, ipak, marr, shape2, shit)

save.image(paste0('mall', '011225.Rdata'))
#save.image("MALL_hex.Rdata")
