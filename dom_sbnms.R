#This script calculates overlap between shearwaters and sightings of HUWH in SBNMS
#uses combination of my code and KP code
#11/8/21

library(adehabitatHR) #can export shape files of polygons
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(tidyverse)
library(readxl)
library(ggpubr)
library(marmap)

#set up stuff for projections and mapping

#proj_NAD83<- CRS("+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")
proj_wgs84 <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
proj_utm <- CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#read in sbnms shapefile
sbnms<-readOGR("./data/sbnms.shp") 
sbnms <- spTransform(sbnms, proj_wgs84)
sbnms_fort <- fortify(sbnms) #for ggplot

#read in 40m contou shapefile
#bath40<-readOGR("./data/bathy40m.shp") 
#bath40 <- spTransform(bath40, proj_wgs84)
#bath40_fort <- fortify(bath40) #for ggplot

#get bathy from marmap so I can have contours
bath_cont <- getNOAA.bathy (lon1 = -65, lon2 = -72.5, lat1 = 40, lat2 = 45, resolution = 1)
bath_cont_fort<- fortify.bathy(bath_cont) #for ggplot


###########################################################

#load in and format whale data

###whale stuff
### load in CCS sightings data
sightings <- read_excel("./data/SBNMS_SIGHTING_LOCATIONS_CCS.xlsx",na=" ") 

#formatting
colnames(sightings) <- c("date","whale_id","lat","lon", "age_class", "repro_state")

#replace . in date character with a - so I can format as date, I don't know why this works
#https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
sightings$date <- gsub('\\.', '-', sightings$date)

#format date as dates
sightings$date <- as.Date(sightings$date,"%Y-%m-%d")

#create year, month and day column
sightings <- sightings %>%
  separate(date, sep="-", into = c("year", "month", "day"))

#separate takes away the full date column, add it back in just in case
sightings <- mutate(sightings,date = paste(year, month, day, sep = "-"))

#format date as dates again
sightings$date <- as.Date(sightings$date,"%Y-%m-%d")

#format month as month
sightings$month <- as.factor(sightings$month)

#filter sightings to july - nov to match bird tag time period
sightings_7to10 <- sightings %>%
  filter(month=="07" | month=="08" | month=="09" | month=="10") %>% droplevels()

#get some sighting sumary information
sightings_7to10 %>% group_by(year) %>% summarise(sightperyear = n())
numsight<- sightings_7to10 %>% group_by(year, month) %>% summarise(sightperyearmonth = n())

###########################################################################################


#2013

#get whale sightings

sight13 <- sightings_7to10 %>% filter(year=='2013')

#get in utm

sight13_xy <- data.frame(sight13$lon, sight13$lat)
coordinates(sight13_xy) <- c("sight13.lon", "sight13.lat")
proj4string(sight13_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight13_utm <- spTransform(sight13_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight13_utm, "SpatialPoints")

#load in files
#telemetry years have different processing bc using combined GRSH-HUWH files needed to calculate overlap of telemetry UDs

#2013

#combined file 
birds_whales13 <- read_csv("./data/huwh-grsh GOM hulls (2013).csv", col_names=TRUE)

birds_whales13 <- birds_whales13 %>% dplyr::select(id,species,dt,utm.x,utm.y,coords.x, coords.y,year,month)
colnames(birds_whales13) <- c("tag", "spp", "datetime","x","y", "lon", "lat", "year","month")

#filter locations to GOM (same KP used): lat < -65, lon 35-49
bw13 <- filter(birds_whales13, lat > 40 & lat < 45) 
bw13 <- filter(bw13, lon < -65)


#filter tracks to just july-october
bw13 <- filter(bw13, month!="11") %>% droplevels()


bw13$tag <- as.factor(bw13$tag)

bw13orig <- bw13

##some summary values

#get trackpts and idivdiuals by year
bw13 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())

# #unique animalsby year
#bw13 %>% group_by(spp,year) %>% summarise(count=n_distinct(tag))

# # #get track pts for each month
#bw13 %>% group_by(spp,year, month) %>% summarise(ptsperyearmonth = n())

# #unique animals by month
#bw13 %>% group_by(spp,year, month) %>% summarise(count=n_distinct(tag))


#get locations ready for UD calc

coordinates(bw13)<- ~x+y		# set the projection on UTM coords
proj4string(bw13)<- proj_utm

#calculae kernel UD - using settings developed by KP
kudl13 <- kernelUD(bw13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                   same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                   boundary = NULL)


#kernel.area(kudl13, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudl13, method = c("HR"), percent = 50, conditional = FALSE)

bw13_ud50 <- getverticeshr(kudl13, percent = 50)

#filter to get just bird ud

bw13_ud50gs <-subset(bw13_ud50, bw13_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight13_in_poly <- over(sight13_utm, bw13_ud50gs)

#find points outside polygon
points13_outside_all_tracks<- sum(!complete.cases(sight13_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap13 <- (nrow(sight13) - points13_outside_all_tracks) / nrow(sight13_in_poly) * 100


#for plotting
#convert to wgs84 for plotting
bw13_ud50_xy <- spTransform(bw13_ud50gs, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
 
# #fortify objects for ggplot

bw13_ud50_fort <-fortify(bw13_ud50_xy)

#need separate data bc of this stuuuuupid code
gs13 <-filter(bw13_ud50_fort, id=="GRSH")

########################################################################################
###months

#2013

#July 2013

#whale data

#subset whale sightings
sight13jul<- subset(sight13, sight13$month == '07')

#get in utm

sight13jul_xy <- data.frame(sight13jul$lon, sight13jul$lat)
coordinates(sight13jul_xy) <- c("sight13jul.lon", "sight13jul.lat")
proj4string(sight13jul_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight13jul_utm <- spTransform(sight13jul_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight13jul_utm, "SpatialPoints")


#birds

jul13<- subset(bw13, bw13$month == 7)

#calculae kernel UD - using settings developed by KP
kudljul13 <- kernelUD(jul13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudljul13, percent = 50, 
 #           unin = c("m"),
 #           unout = c("km2"))

#kerneloverlaphr(kudljul13, method = c("HR"), percent = 50, conditional = FALSE)

jul13_ud50 <- getverticeshr(kudljul13, percent = 50)

#filter to get just bird ud

jul13_ud50gs <-subset(jul13_ud50, jul13_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight13jul_in_poly <- over(sight13jul_utm, jul13_ud50gs)

#find points outside polygon
points13jul_outside_all_tracks<- sum(!complete.cases(sight13jul_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap13jul <- (nrow(sight13jul) - points13jul_outside_all_tracks) / nrow(sight13jul_in_poly) * 100


#for plotting
#convert to wgs84 for plotting
jul13_ud50_xy <- spTransform(jul13_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
jul13_ud50_fort <-fortify(jul13_ud50_xy)

#need separate data bc of this stuuuuupid code
gsjul13 <-filter(jul13_ud50_fort, id=="GRSH")


###aug 2013

#whale data

#subset whale sightings
sight13aug<- subset(sight13, sight13$month == '08')

#get in utm

sight13aug_xy <- data.frame(sight13aug$lon, sight13aug$lat)
coordinates(sight13aug_xy) <- c("sight13aug.lon", "sight13aug.lat")
proj4string(sight13aug_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight13aug_utm <- spTransform(sight13aug_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight13aug_utm, "SpatialPoints")


aug13<- subset(bw13, bw13$month == 8)

#calculae kernel UD - using settings developed by KP
kudlaug13 <- kernelUD(aug13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudlaug13, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudlaug13, method = c("HR"), percent = 50, conditional = FALSE)

aug13_ud50 <- getverticeshr(kudlaug13, percent = 50)

#filter to get just bird ud

aug13_ud50gs <-subset(aug13_ud50, aug13_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight13aug_in_poly <- over(sight13aug_utm, aug13_ud50gs)

#find points outside polygon
points13aug_outside_all_tracks<- sum(!complete.cases(sight13aug_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap13aug <- (nrow(sight13aug) - points13aug_outside_all_tracks) / nrow(sight13aug_in_poly) * 100

#for plotting
#convert to wgs84 for plotting
aug13_ud50_xy <- spTransform(aug13_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
aug13_ud50_fort <-fortify(aug13_ud50_xy)

#need separate data bc of this stuuuuupid code
gsaug13 <-filter(aug13_ud50_fort, id=="GRSH")



###sep 2013

#whale data

#subset whale sightings
sight13sep<- subset(sight13, sight13$month == '09')

#get in utm

sight13sep_xy <- data.frame(sight13sep$lon, sight13sep$lat)
coordinates(sight13sep_xy) <- c("sight13sep.lon", "sight13sep.lat")
proj4string(sight13sep_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight13sep_utm <- spTransform(sight13sep_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight13sep_utm, "SpatialPoints")


sep13<- subset(bw13, bw13$month == 9)

#calculae kernel UD - using settings developed by KP
kudlsep13 <- kernelUD(sep13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudlsep13, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudlsep13, method = c("HR"), percent = 50, conditional = FALSE)

sep13_ud50 <- getverticeshr(kudlsep13, percent = 50)

#filter to get just bird ud

sep13_ud50gs <-subset(sep13_ud50, sep13_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight13sep_in_poly <- over(sight13sep_utm, sep13_ud50gs)

#find points outside polygon
points13sep_outside_all_tracks<- sum(!complete.cases(sight13sep_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap13sep <- (nrow(sight13sep) - points13sep_outside_all_tracks) / nrow(sight13sep_in_poly) * 100


#for plotting
#convert to wgs84 for plotting
sep13_ud50_xy <- spTransform(sep13_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
sep13_ud50_fort <-fortify(sep13_ud50_xy)

#need separate data bc of this stuuuuupid code
gssep13 <-filter(sep13_ud50_fort, id=="GRSH")



###oct 2013

#whale data

#subset whale sightings
sight13oct<- subset(sight13, sight13$month == '10')

#get in utm

sight13oct_xy <- data.frame(sight13oct$lon, sight13oct$lat)
coordinates(sight13oct_xy) <- c("sight13oct.lon", "sight13oct.lat")
proj4string(sight13oct_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight13oct_utm <- spTransform(sight13oct_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight13oct_utm, "SpatialPoints")


oct13<- subset(bw13, bw13$month == 10)

#calculae kernel UD - using settings developed by KP
kudloct13 <- kernelUD(oct13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudloct13, percent = 50, 
#            unin = c("m"),
 #           unout = c("km2"))

#kerneloverlaphr(kudloct13, method = c("HR"), percent = 50, conditional = FALSE)

oct13_ud50 <- getverticeshr(kudloct13, percent = 50)

#filter to get just bird ud

oct13_ud50gs <-subset(oct13_ud50, oct13_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight13oct_in_poly <- over(sight13oct_utm, oct13_ud50gs)

#find points outside polygon
points13oct_outside_all_tracks<- sum(!complete.cases(sight13oct_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap13oct <- (nrow(sight13oct) - points13oct_outside_all_tracks) / nrow(sight13oct_in_poly) * 100


#for plotting
#convert to wgs84 for plotting
oct13_ud50_xy <- spTransform(oct13_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
oct13_ud50_fort <-fortify(oct13_ud50_xy)

#need separate data bc of this stuuuuupid code
gsoct13 <-filter(oct13_ud50_fort, id=="GRSH")


###2015

#get whale sightings

sight15 <- sightings_7to10 %>% filter(year=='2015')

#get in utm

sight15_xy <- data.frame(sight15$lon, sight15$lat)
coordinates(sight15_xy) <- c("sight15.lon", "sight15.lat")
proj4string(sight15_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight15_utm <- spTransform(sight15_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight15_utm, "SpatialPoints")


#bird stuff
#combined file 
birds_whales15 <- read_csv("./data/huwh-grsh GOM hulls (2015).csv", col_names=TRUE)


birds_whales15 <- birds_whales15 %>% dplyr::select(id,species,dt,utm.x,utm.y,coords.x, coords.y,year,month)
colnames(birds_whales15) <- c("tag", "spp", "datetime","x","y", "lon", "lat", "year","month")


#lat and lon are reversed for grsh, need to swap cols

#pull out and save locations
grsh <- birds_whales15 %>% filter(spp =='GRSH')
grsh_lat<- grsh$lat
grsh_lon <- grsh$lon

#find indicies in df = GRSH
range(which(birds_whales15$spp == 'GRSH'))

birds_whales15$lat[14425:23598] <- grsh_lon

birds_whales15$lon[14425:23598] <- grsh_lat

#filter locations to GOM (same KP used): lat < -65, lon 35-49
bw15 <- filter(birds_whales15, lat > 40 & lat < 45) 
bw15 <- filter(bw15, lon < -65)


#filter tracks to just july-october
bw15 <- filter(bw15, month!="11") %>% droplevels()


bw15$tag <- as.factor(bw15$tag)

bw15orig <- bw15


##some summary values

#get trackpts and idivdiuals by year
#bw15 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animalsby year
#bw15 %>% group_by(spp,year) %>% summarise(count=n_distinct(tag))
# 
# # #get track pts for each month
#bw15 %>% group_by(spp,year, month) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animals by month
#bw15 %>% group_by(spp,year, month) %>% summarise(count=n_distinct(tag))


#get locations ready for UD calc

coordinates(bw15)<- ~x+y		# set the projection on UTM coords
proj4string(bw15)<- proj_utm

#calculae kernel UD - using settings developed by KP
kudl15 <- kernelUD(bw15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                   same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                   boundary = NULL)


#kernel.area(kudl15, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudl15, method = c("HR"), percent = 50, conditional = FALSE)

bw15_ud50 <- getverticeshr(kudl15, percent = 50)


#filter to get just bird ud

bw15_ud50gs <-subset(bw15_ud50, bw15_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight15_in_poly <- over(sight15_utm, bw15_ud50gs)

#find points outside polygon
points15_outside_all_tracks<- sum(!complete.cases(sight15_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap15 <- (nrow(sight15) - points15_outside_all_tracks) / nrow(sight15_in_poly) * 100

percent_overlap15

#for plotting
#convert to wgs84 for plotting
bw15_ud50_xy <- spTransform(bw15_ud50gs, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

# #fortify objects for ggplot

bw15_ud50_fort <-fortify(bw15_ud50_xy)

#need separate data bc of this stuuuuupid code
gs15 <-filter(bw15_ud50_fort, id=="GRSH")

#for plotting
#convert to wgs84 for plotting
bw15_ud50_xy <- spTransform(bw15_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
bw15_ud50_fort <-fortify(bw15_ud50_xy)

#need separate data bc of this stuuuuupid code
gs15 <-filter(bw15_ud50_fort, id=="GRSH")
hw15 <-filter(bw15_ud50_fort, id=="HUWH")

###months
#2015

#July 2015

#whale data

#subset whale sightings
sight15jul<- subset(sight15, sight15$month == '07')

#get in utm

sight15jul_xy <- data.frame(sight15jul$lon, sight15jul$lat)
coordinates(sight15jul_xy) <- c("sight15jul.lon", "sight15jul.lat")
proj4string(sight15jul_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight15jul_utm <- spTransform(sight15jul_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight15jul_utm, "SpatialPoints")


#birds

jul15<- subset(bw15, bw15$month == 7)

#calculae kernel UD - using settings developed by KP
kudljul15 <- kernelUD(jul15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudljul15, percent = 50, 
#           unin = c("m"),
#           unout = c("km2"))

#kerneloverlaphr(kudljul15, method = c("HR"), percent = 50, conditional = FALSE)

jul15_ud50 <- getverticeshr(kudljul15, percent = 50)

#filter to get just bird ud

jul15_ud50gs <-subset(jul15_ud50, jul15_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight15jul_in_poly <- over(sight15jul_utm, jul15_ud50gs)

#find points outside polygon
points15jul_outside_all_tracks<- sum(!complete.cases(sight15jul_in_poly))


#find percentage of sightings inside the all_tracks ud50
percent_overlap15jul <- (nrow(sight15jul) - points15jul_outside_all_tracks) / nrow(sight15jul_in_poly) * 100

percent_overlap15jul

#for plotting
#convert to wgs84 for plotting
jul15_ud50_xy <- spTransform(jul15_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
jul15_ud50_fort <-fortify(jul15_ud50_xy)

#need separate data bc of this stuuuuupid code
gsjul15 <-filter(jul15_ud50_fort, id=="GRSH")


###aug 2015

#whale data

#subset whale sightings
sight15aug<- subset(sight15, sight15$month == '08')

#get in utm

sight15aug_xy <- data.frame(sight15aug$lon, sight15aug$lat)
coordinates(sight15aug_xy) <- c("sight15aug.lon", "sight15aug.lat")
proj4string(sight15aug_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight15aug_utm <- spTransform(sight15aug_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight15aug_utm, "SpatialPoints")


aug15<- subset(bw15, bw15$month == 8)

#calculae kernel UD - using settings developed by KP
kudlaug15 <- kernelUD(aug15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudlaug15, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudlaug15, method = c("HR"), percent = 50, conditional = FALSE)

aug15_ud50 <- getverticeshr(kudlaug15, percent = 50)

#filter to get just bird ud

aug15_ud50gs <-subset(aug15_ud50, aug15_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight15aug_in_poly <- over(sight15aug_utm, aug15_ud50gs)

#find points outside polygon
points15aug_outside_all_tracks<- sum(!complete.cases(sight15aug_in_poly))

points15aug_outside_all_tracks


#find percentage of sightings inside the all_tracks ud50
percent_overlap15aug <- (nrow(sight15aug) - points15aug_outside_all_tracks) / nrow(sight15aug_in_poly) * 100

percent_overlap15aug

#for plotting
#convert to wgs84 for plotting
aug15_ud50_xy <- spTransform(aug15_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
aug15_ud50_fort <-fortify(aug15_ud50_xy)

#need separate data bc of this stuuuuupid code
gsaug15 <-filter(aug15_ud50_fort, id=="GRSH")



###sep 2015

#whale data

#subset whale sightings
sight15sep<- subset(sight15, sight15$month == '09')

#get in utm

sight15sep_xy <- data.frame(sight15sep$lon, sight15sep$lat)
coordinates(sight15sep_xy) <- c("sight15sep.lon", "sight15sep.lat")
proj4string(sight15sep_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight15sep_utm <- spTransform(sight15sep_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight15sep_utm, "SpatialPoints")


sep15<- subset(bw15, bw15$month == 9)

#calculae kernel UD - using settings developed by KP
kudlsep15 <- kernelUD(sep15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudlsep15, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudlsep15, method = c("HR"), percent = 50, conditional = FALSE)

sep15_ud50 <- getverticeshr(kudlsep15, percent = 50)

#filter to get just bird ud

sep15_ud50gs <-subset(sep15_ud50, sep15_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight15sep_in_poly <- over(sight15sep_utm, sep15_ud50gs)

#find points outside polygon
points15sep_outside_all_tracks<- sum(!complete.cases(sight15sep_in_poly))

points15sep_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap15sep <- (nrow(sight15sep) - points15sep_outside_all_tracks) / nrow(sight15sep_in_poly) * 100

percent_overlap15sep

#for plotting
#convert to wgs84 for plotting
sep15_ud50_xy <- spTransform(sep15_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
sep15_ud50_fort <-fortify(sep15_ud50_xy)

#need separate data bc of this stuuuuupid code
gssep15 <-filter(sep15_ud50_fort, id=="GRSH")



###oct 2015

#whale data

#subset whale sightings
sight15oct<- subset(sight15, sight15$month == '10')

#get in utm

sight15oct_xy <- data.frame(sight15oct$lon, sight15oct$lat)
coordinates(sight15oct_xy) <- c("sight15oct.lon", "sight15oct.lat")
proj4string(sight15oct_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight15oct_utm <- spTransform(sight15oct_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight15oct_utm, "SpatialPoints")


oct15<- subset(bw15, bw15$month == 10)

#calculae kernel UD - using settings developed by KP
kudloct15 <- kernelUD(oct15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudloct15, percent = 50, 
#            unin = c("m"),
#           unout = c("km2"))

#kerneloverlaphr(kudloct15, method = c("HR"), percent = 50, conditional = FALSE)

oct15_ud50 <- getverticeshr(kudloct15, percent = 50)

#filter to get just bird ud

oct15_ud50gs <-subset(oct15_ud50, oct15_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight15oct_in_poly <- over(sight15oct_utm, oct15_ud50gs)

#find points outside polygon
points15oct_outside_all_tracks<- sum(!complete.cases(sight15oct_in_poly))

points15oct_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap15oct <- (nrow(sight15oct) - points15oct_outside_all_tracks) / nrow(sight15oct_in_poly) * 100

percent_overlap15oct

#for plotting
#convert to wgs84 for plotting
oct15_ud50_xy <- spTransform(oct15_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
oct15_ud50_fort <-fortify(oct15_ud50_xy)

#need separate data bc of this stuuuuupid code
gsoct15 <-filter(oct15_ud50_fort, id=="GRSH")


###2018

#whale data

#get whale sightings

sight18 <- sightings_7to10 %>% filter(year=='2018')

#get in utm

sight18_xy <- data.frame(sight18$lon, sight18$lat)
coordinates(sight18_xy) <- c("sight18.lon", "sight18.lat")
proj4string(sight18_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight18_utm <- spTransform(sight18_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight18_utm, "SpatialPoints")

#combined file 
birds_whales18 <- read_csv("./data/huwh-grsh GOM hulls (2018).csv", col_names=TRUE)


birds_whales18 <- birds_whales18 %>% dplyr::select(id,species,dt,utm.x,utm.y,coords.x, coords.y,year,month)
colnames(birds_whales18) <- c("tag", "spp", "datetime","x","y", "lon", "lat", "year","month")


#lat and lon are reversed for grsh, need to swap cols

#pull out and save locations
grsh <- birds_whales18 %>% filter(spp =='GRSH')
grsh_lat<- grsh$lat
grsh_lon <- grsh$lon

#find indicies in df = GRSH
range(which(birds_whales18$spp == 'GRSH'))

birds_whales18$lat[5420:13395] <- grsh_lon

birds_whales18$lon[5420:13395] <- grsh_lat

#filter locations to GOM (same KP used): lat < -65, lon 35-49
bw18 <- filter(birds_whales18, lat > 40 & lat < 45) 
bw18 <- filter(bw18, lon < -65)


#filter tracks to just july-october
bw18 <- filter(bw18, month!="11") %>% droplevels()


bw18$tag <- as.factor(bw18$tag)

bw18orig <- bw18

##some summary values

#get trackpts and idivdiuals by year
#bw18 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animalsby year
#bw18 %>% group_by(spp,year) %>% summarise(count=n_distinct(tag))
# 
# # #get track pts for each month
#bw18 %>% group_by(spp,year, month) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animals by month
#bw18 %>% group_by(spp,year, month) %>% summarise(count=n_distinct(tag))


#get locations ready for UD calc

coordinates(bw18)<- ~x+y		# set the projection on UTM coords
proj4string(bw18)<- proj_utm

#calculae kernel UD - using settings developed by KP
kudl18 <- kernelUD(bw18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                   same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                   boundary = NULL)


#kernel.area(kudl18, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudl18, method = c("HR"), percent = 50, conditional = FALSE)

bw18_ud50 <- getverticeshr(kudl18, percent = 50)

#filter to get just bird ud

bw18_ud50gs <-subset(bw18_ud50, bw18_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight18_in_poly <- over(sight18_utm, bw18_ud50gs)

#find points outside polygon
points18_outside_all_tracks<- sum(!complete.cases(sight18_in_poly))

points18_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap18 <- (nrow(sight18) - points18_outside_all_tracks) / nrow(sight18_in_poly) * 100

percent_overlap18


#for plotting
#convert to wgs84 for plotting
bw18_ud50_xy <- spTransform(bw18_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
bw18_ud50_fort <-fortify(bw18_ud50_xy)

#need separate data bc of this stuuuuupid code
gs18 <-filter(bw18_ud50_fort, id=="GRSH")



###months
#2018

#July 2018

#whale data

#subset whale sightings
sight18jul<- subset(sight18, sight18$month == '07')

#get in utm

sight18jul_xy <- data.frame(sight18jul$lon, sight18jul$lat)
coordinates(sight18jul_xy) <- c("sight18jul.lon", "sight18jul.lat")
proj4string(sight18jul_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight18jul_utm <- spTransform(sight18jul_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight18jul_utm, "SpatialPoints")


#birds

jul18<- subset(bw18, bw18$month == 7)

#calculae kernel UD - using settings developed by KP
kudljul18 <- kernelUD(jul18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudljul18, percent = 50, 
#           unin = c("m"),
#           unout = c("km2"))

#kerneloverlaphr(kudljul18, method = c("HR"), percent = 50, conditional = FALSE)

jul18_ud50 <- getverticeshr(kudljul18, percent = 50)

#filter to get just bird ud

jul18_ud50gs <-subset(jul18_ud50, jul18_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight18jul_in_poly <- over(sight18jul_utm, jul18_ud50gs)

#find points outside polygon
points18jul_outside_all_tracks<- sum(!complete.cases(sight18jul_in_poly))

points18jul_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap18jul <- (nrow(sight18jul) - points18jul_outside_all_tracks) / nrow(sight18jul_in_poly) * 100

percent_overlap18jul

#for plotting
#convert to wgs84 for plotting
jul18_ud50_xy <- spTransform(jul18_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
jul18_ud50_fort <-fortify(jul18_ud50_xy)

#need separate data bc of this stuuuuupid code
gsjul18 <-filter(jul18_ud50_fort, id=="GRSH")


###aug 2018

#whale data

#subset whale sightings
sight18aug<- subset(sight18, sight18$month == '08')

#get in utm

sight18aug_xy <- data.frame(sight18aug$lon, sight18aug$lat)
coordinates(sight18aug_xy) <- c("sight18aug.lon", "sight18aug.lat")
proj4string(sight18aug_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight18aug_utm <- spTransform(sight18aug_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight18aug_utm, "SpatialPoints")


aug18<- subset(bw18, bw18$month == 8)

#calculae kernel UD - using settings developed by KP
kudlaug18 <- kernelUD(aug18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudlaug18, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudlaug18, method = c("HR"), percent = 50, conditional = FALSE)

aug18_ud50 <- getverticeshr(kudlaug18, percent = 50)

#filter to get just bird ud

aug18_ud50gs <-subset(aug18_ud50, aug18_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight18aug_in_poly <- over(sight18aug_utm, aug18_ud50gs)

#find points outside polygon
points18aug_outside_all_tracks<- sum(!complete.cases(sight18aug_in_poly))

points18aug_outside_all_tracks


#find percentage of sightings inside the all_tracks ud50
percent_overlap18aug <- (nrow(sight18aug) - points18aug_outside_all_tracks) / nrow(sight18aug_in_poly) * 100

percent_overlap18aug

#for plotting
#convert to wgs84 for plotting
aug18_ud50_xy <- spTransform(aug18_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
aug18_ud50_fort <-fortify(aug18_ud50_xy)

#need separate data bc of this stuuuuupid code
gsaug18 <-filter(aug18_ud50_fort, id=="GRSH")



###sep 2018

#whale data

#subset whale sightings
sight18sep<- subset(sight18, sight18$month == '09')

#get in utm

sight18sep_xy <- data.frame(sight18sep$lon, sight18sep$lat)
coordinates(sight18sep_xy) <- c("sight18sep.lon", "sight18sep.lat")
proj4string(sight18sep_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight18sep_utm <- spTransform(sight18sep_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight18sep_utm, "SpatialPoints")


sep18<- subset(bw18, bw18$month == 9)

#calculae kernel UD - using settings developed by KP
kudlsep18 <- kernelUD(sep18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudlsep18, percent = 50, 
#            unin = c("m"),
#            unout = c("km2"))

#kerneloverlaphr(kudlsep18, method = c("HR"), percent = 50, conditional = FALSE)

sep18_ud50 <- getverticeshr(kudlsep18, percent = 50)

#filter to get just bird ud

sep18_ud50gs <-subset(sep18_ud50, sep18_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight18sep_in_poly <- over(sight18sep_utm, sep18_ud50gs)

#find points outside polygon
points18sep_outside_all_tracks<- sum(!complete.cases(sight18sep_in_poly))

points18sep_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap18sep <- (nrow(sight18sep) - points18sep_outside_all_tracks) / nrow(sight18sep_in_poly) * 100

percent_overlap18sep

#for plotting
#convert to wgs84 for plotting
sep18_ud50_xy <- spTransform(sep18_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
sep18_ud50_fort <-fortify(sep18_ud50_xy)

#need separate data bc of this stuuuuupid code
gssep18 <-filter(sep18_ud50_fort, id=="GRSH")



###oct 2018

#whale data

#subset whale sightings
sight18oct<- subset(sight18, sight18$month == '10')

#get in utm

sight18oct_xy <- data.frame(sight18oct$lon, sight18oct$lat)
coordinates(sight18oct_xy) <- c("sight18oct.lon", "sight18oct.lat")
proj4string(sight18oct_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sight18oct_utm <- spTransform(sight18oct_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sight18oct_utm, "SpatialPoints")


oct18<- subset(bw18, bw18$month == 10)

#calculae kernel UD - using settings developed by KP
kudloct18 <- kernelUD(oct18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


#kernel.area(kudloct18, percent = 50, 
#            unin = c("m"),
#           unout = c("km2"))

#kerneloverlaphr(kudloct18, method = c("HR"), percent = 50, conditional = FALSE)

oct18_ud50 <- getverticeshr(kudloct18, percent = 50)

#filter to get just bird ud

oct18_ud50gs <-subset(oct18_ud50, oct18_ud50$id=="GRSH")


###get overlap

#figure out how many sightings fall within the all_tracks ud50
sight18oct_in_poly <- over(sight18oct_utm, oct18_ud50gs)

#find points outside polygon
points18oct_outside_all_tracks<- sum(!complete.cases(sight18oct_in_poly))

points18oct_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap18oct <- (nrow(sight18oct) - points18oct_outside_all_tracks) / nrow(sight18oct_in_poly) * 100

percent_overlap18oct

#for plotting
#convert to wgs84 for plotting
oct18_ud50_xy <- spTransform(oct18_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
oct18_ud50_fort <-fortify(oct18_ud50_xy)

#need separate data bc of this stuuuuupid code
gsoct18 <-filter(oct18_ud50_fort, id=="GRSH")


###2014 to 2017

tracks14 <- read_excel("./data/tlocoh_lhs_locs (2014).xlsx")
tracks16 <- read_excel("./data/tlocoh_lhs_locs (2016).xlsx")
tracks17 <- read_excel("./data/tlocoh_lhs_locs (2017).xlsx")

tracks14 <- tracks14 %>% dplyr::select(id,dt,x,y,coords.x, coords.y,year,month)
tracks16 <- tracks16 %>% dplyr::select(id,dt,x,y,coords.x, coords.y,year,month)
tracks17 <- tracks17 %>% dplyr::select(id,dt,x,y,coords.x, coords.y,year,month)

colnames(tracks14) <- c("tag", "datetime","x","y", "lon", "lat", "year","month")
colnames(tracks16) <- c("tag", "datetime","x","y", "lon", "lat", "year","month")
colnames(tracks17) <- c("tag", "datetime","x","y", "lon", "lat", "year","month")

#need to remove sooty shearwater Brady
tracks17 <- tracks17 %>% filter(tag !='Brady')

#combine all data
all_tracks <- rbind(tracks14,tracks16,tracks17)

#format datetime as GMT
attr(all_tracks$datetime, "tzone") <- 'GMT'

#make a local datetime (raw data is UTC)
all_tracks$local_datetime <- all_tracks$datetime
attr(all_tracks$local_datetime, "tzone") <- 'EST'

#filter locations to GOM (same KP used): lat < -65, lon 35-49
all_tracks_gom <- filter(all_tracks, lat > 40 & lat < 45) 
all_tracks_gom <- filter(all_tracks_gom, lon < -65)

#filter tracks to just july-october
all_tracks_gom <- filter(all_tracks_gom, month!="11") %>% droplevels()

#split by year and make a list and get in utm

#whale sightings
#first just use 2014, 2016-2017
sightings14to17 <- sightings_7to10 %>% filter(year=='2014' | year=='2016' | year=='2017') %>% droplevels()
dat_list <- split(sightings14to17, sightings14to17$year) 
list_xy<- lapply(dat_list, function(x){
  sightings_xy <- data.frame(x$lon, x$lat)
  coordinates(sightings_xy) <- c("x.lon", "x.lat")
  proj4string(sightings_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  sightings_utm <- spTransform(sightings_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(sightings_utm, "SpatialPoints")
})

bird_list <- split(all_tracks_gom, all_tracks_gom$year) 
bird_list_xy<- lapply(bird_list, function(y){
  bird_xy <- data.frame(y$lon, y$lat)
  coordinates(bird_xy) <- c("y.lon", "y.lat")
  proj4string(bird_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  bird_utm <- spTransform(bird_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(bird_utm, "SpatialPoints")
})

all_tracks_kd <- lapply(bird_list_xy, function(y){
  kernelUD(y ,h="href", grid = 80, kern = c("bivnorm"), 
           same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
           boundary = NULL)#estimate 90 and 50% UD for july2014 GOM data
})

#find the area of each 50%UD
areas <- lapply(all_tracks_kd, function(y){
  kernel.area(y, percent = 50,
              unin = c("m"),
              unout = c("km2"))
})

areas

all_tracks_kd <- lapply(bird_list_xy, function(y){
  kud <- kernelUD(y ,h="href", grid = 80, kern = c("bivnorm"), 
                  same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                  boundary = NULL)#estimate 90 and 50% UD for july2014 GOM data
  getverticeshr(kud,percent=50)#save the ud as a spatialpolygonsdataframe
})


sightings_outside <- mapply(function(x,y){
  sightings_in_poly <- over(x,y)
  sightings_outside_poly <- sum(!complete.cases(sightings_in_poly))
},
list_xy, all_tracks_kd) 

#make overlap a dataframe -this is number of sightings outside the poly
sightings_outside <- data.frame(sightings_outside)
colnames(sightings_outside) <- c("points_outside")

#make year a column in df
sightings_outside <- sightings_outside %>%
  rownames_to_column(var="year")

#split into list by year
sightings_outside <- split(sightings_outside,sightings_outside$year)

percent_overlap <- mapply(function(x,y){
  percent_overlap <- (nrow(x) - y$points_outside) / nrow(x) * 100
},
dat_list,sightings_outside)  

percent_overlap #get percent overlap for each year - percentage of  huwh sightings inside grsh 50ud

#get the number of trackpoints and tagged birds used for analysis in each year
all_tracks_gom %>% group_by(year) %>% summarise(trackptsperyear = n(), indivi = n_distinct(tag))

#for plotting
all_tracks_kd_xy <- lapply(all_tracks_kd, function(x){
  spTransform(x, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
})

#fortfy objects for ggplot
all_tracks_kd_xy_fort <- lapply(all_tracks_kd_xy, function(x){
  fortify(x)
})

###Months 14 to 17

#remove July from 2017 which has no trackpoints
sightings14to17 <- sightings14to17 %>%
  filter(!(year=='2017' & month=='07')) %>%
  droplevels()

#split dataframe by year and by month
dat_list_month <- group_split(sightings14to17, year, month)

list_xy_month <- lapply(dat_list_month, function(x){
  sightings_xy_month <- data.frame(x$lon, x$lat)
  coordinates(sightings_xy_month) <- c("x.lon", "x.lat")
  proj4string(sightings_xy_month) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  sightings_utm_month <- spTransform(sightings_xy_month, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(sightings_utm_month, "SpatialPoints")
})

bird_list_month <- group_split(all_tracks_gom, year, month)

bird_list_xy_month<- lapply(bird_list_month, function(y){
  bird_xy_month <- data.frame(y$lon, y$lat)
  coordinates(bird_xy_month) <- c("y.lon", "y.lat")
  proj4string(bird_xy_month) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  bird_utm_month <- spTransform(bird_xy_month, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(bird_utm_month, "SpatialPoints")
})

#take the new bird list and do all the kernel density stuff

all_tracks_kd <- lapply(bird_list_xy_month, function(y){
  kernelUD(y ,h="href", grid = 80, kern = c("bivnorm"), 
           same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
           boundary = NULL)#estimate 90 and 50% UD for july2014 GOM data
})

#find the area of each 50%UD
areas <- lapply(all_tracks_kd, function(y){ #figure out how to get month labels on here
  kernel.area(y, percent = 50,
              unin = c("m"),
              unout = c("km2"))
})

#save to csv
write.csv(areas, "areas_months14to17.csv") 

all_tracks_kd_month <- lapply(bird_list_xy_month, function(y){
  kud_month <- kernelUD(y ,h="href", grid = 80, kern = c("bivnorm"), 
                        same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                        boundary = NULL) #estimate 90 and 50% UD for july2014 GOM data
  getverticeshr(kud_month,percent=50)#save the ud as a spatialpolygonsdataframe
})

#now I have a bird list with the hulls and a whale list in utm
#find sightings in the hulls

sightings_outside_month <- mapply(function(x,y){
  sightings_in_poly_month <- over(x,y)
  sightings_outside_poly_month <- sum(!complete.cases(sightings_in_poly_month))
},
list_xy_month, all_tracks_kd_month) 

#make overlap a dataframe -this is number of sightings outside the poly
sightings_outside_month <- data.frame(sightings_outside_month)
colnames(sightings_outside_month) <- c("points_outside")
sightings_outside_month <- split(sightings_outside_month,seq(nrow(sightings_outside_month))) #had to use seq(nrow) bc for some reason just split was dropping rows with 0 values

percent_overlap_month <- mapply(function(x,y){
  percent_overlap_month <- (nrow(x) - y$points_outside) / nrow(x) * 100
},
dat_list_month,sightings_outside_month)  

percent_overlap_month #get percent overlap for each year - percentage of  huwh sightings inside grsh 50ud

#get the number of trackpoints used for analysis in each year

trackptsyrmonth<- all_tracks_gom %>% group_by(year, month) %>% summarise(trackptsperyearmonth = n(), tags=n_distinct(tag))

#for ploting
#convert to wgs84 for plotting
all_tracks_kd_xy_month <- lapply(all_tracks_kd_month, function(x){
  spTransform(x, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
})

#fortfy objects for ggplot
all_tracks_kd_xy_month_fort <- lapply(all_tracks_kd_xy_month, function(x){
  fortify(x)
})




#####################################################################################################

###Plotting

#years

world_map <- map_data("world")

#2013

p13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight13, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gs13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=gs13,location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.1,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = gs13, location = "bottomright", scale=0.4, symbol=1, anchor=c(x=-65, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  


#2014

p14 <-  ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list[[1]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_fort[[1]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 



#2015

p15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight15, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gs15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 


#2016

p16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list[[2]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_fort[[2]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "95%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 


# 2017

p17 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list[[3]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_fort[[3]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "91%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))   


# 2018

p18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight18, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gs18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 

#make multipanel plot of annual overlap
#amazing reference / help for all you can do in ggarrange 
#https://www.r-bloggers.com/2017/07/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

all_years <- ggarrange(p13 + rremove("ylab") + rremove("xlab"),
                       p14+ rremove("ylab") + rremove("xlab"),
                       p15+ rremove("ylab") + rremove("xlab"),
                       p16+ rremove("ylab") + rremove("xlab"),
                       p17+ rremove("ylab") + rremove("xlab"),
                       p18 + rremove("ylab") + rremove("xlab"),
                       labels = c("A", "B", "C", "D", "E", "F"),
                       ncol=2,nrow=3) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(all_years, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/all_years13to18_sbnms.tiff",dpi=300,width=6,height=8.5,units="in")

##########################################################################

#months

#2013

#july

j13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight13jul, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsjul13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=gsjul13,location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = gsjul13, location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="July 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


##August 2013


a13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight13aug, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsaug13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


##Sept 2013


s13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight13sep, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gssep13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="September 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

###Plot Oct 2013

o13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight13oct, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsoct13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="October 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


months13 <- ggarrange(j13 + rremove("ylab") + rremove("xlab"),
                      a13+ rremove("ylab") + rremove("xlab"),
                      s13+ rremove("ylab") + rremove("xlab"),
                      o13+ rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months13, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2013_monthly_sbnms.tiff",dpi=300,width=6,height=6,units="in")

###2015

###Plot July 2015

j15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight15jul, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsjul15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=gsjul15,location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = gsjul15, location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="July 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Aug 2015

a15 <- ggplot() +
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight15aug, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsaug15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "95%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

###Plot Sept 2015

s15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight15sep, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gssep15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Oct 2015

o15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight15oct, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsoct15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="October 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months15 <- ggarrange(j15 + rremove("ylab") + rremove("xlab"),
                      a15+ rremove("ylab") + rremove("xlab"),
                      s15+ rremove("ylab") + rremove("xlab"),
                      o15+ rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months15, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2015_monthly_sbnms.tiff",dpi=300,width=6,height=6,units="in")


#2018
#July 2018


j18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight18jul, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsjul18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=gsjul18,location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = gsjul18, location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title=" July 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot August 2018
#quick plot check to confirm %overlap

a18<- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight18aug, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsaug18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Sept 2018
#quick plot check to confirm %overlap

s18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight18sep, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gssep18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Oct 2018

o18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=sight18oct, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = gsoct18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="October 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months18 <- ggarrange(j18 + rremove("ylab") + rremove("xlab"),
                      a18 + rremove("ylab") + rremove("xlab"),
                      s18 + rremove("ylab") + rremove("xlab"),
                      o18 + rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months18, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2018_monthly_sbnms.tiff",dpi=300,width=6,height=6,units="in")

###################################################################################



###Plotting months
#2014-2017
###Plot July 2014

world_map <- map_data("world")
j14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[1]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[1]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=all_tracks_kd_xy_month_fort[[1]],location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = all_tracks_kd_xy_month_fort[[1]], location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="July 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



###Plot August 2014
#quick plot check to confirm %overlap

a14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[2]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[2]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Sept 2014
#quick plot check to confirm %overlap

s14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[3]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[3]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "97%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Oct 2014

o14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[4]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[4]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="October 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months14 <- ggarrange(j14 + rremove("ylab") + rremove("xlab"),
                      a14+ rremove("ylab") + rremove("xlab"),
                      s14+ rremove("ylab") + rremove("xlab"),
                      o14+ rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months14, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2014_monthly_sbnms.tiff",dpi=300,width=6,height=6,units="in")


j16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[5]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
#ud50 for all tracks combined in GOM
geom_polygon(data = all_tracks_kd_xy_month_fort[[5]],
             aes(x = long, y = lat, group = group),inherit.aes = FALSE,
             fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=all_tracks_kd_xy_month_fort[[5]],location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = all_tracks_kd_xy_month_fort[[5]], location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="July 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "90%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

###Plot Aug 2016

a16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[6]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[6]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "91%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Sept 2016

s16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[7]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[7]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot oct 2016

o16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[8]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[8]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="October 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "29%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months16 <- ggarrange(j16 + rremove("ylab") + rremove("xlab"),
                      a16+ rremove("ylab") + rremove("xlab"),
                      s16+ rremove("ylab") + rremove("xlab"),
                      o16+ rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months16, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2016_monthly_sbnms.tiff",dpi=300,width=6,height=6,units="in")


###Plot for July 2017 - just sightings
#refilter to get sightings
# sight717 <- sightings %>% filter(year=='2017' & month=='07')
# 
# j17 <- ggplot() + 
#   #   #add in sightings
#   geom_point(data=sight717, aes(x=lon, y=lat),inherit.aes = FALSE,
#              size=0.1, show.legend=FALSE) +
#   
#   #land
#   geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
#   coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
#   #bathym40m
#   geom_path(data = bath40_fort,
#             aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'darkgray', size = .2) +
#   #sbnms
#   geom_path(data = sbnms_fort,
#             aes(x = long, y = lat, group = group),inherit.aes = FALSE,
#             color = 'black', size = .2) +
#   
#   labs(x="Longitude", y="Latitude",title="July 2017") +
#   annotate("text",x=-70.3,y=44.6,size=5, label = "No tagged birds",color="white") +
#   theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
#                                  plot.title=element_text(size=12))


###Plot Aug 2017

a17 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[9]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[9]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=all_tracks_kd_xy_month_fort[[8]],location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = all_tracks_kd_xy_month_fort[[8]], location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="August 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Sept 2017

s17 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[10]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[10]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "70%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Oct 2017
o17 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_month[[11]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_month_fort[[11]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="October 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "74%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months17 <- ggarrange(
                      a17+ rremove("ylab") + rremove("xlab"),
                      s17+ rremove("ylab") + rremove("xlab"),
                      o17+ rremove("ylab") + rremove("xlab"),
                
                      labels = c("A", "B", "C"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months17, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2017_monthly_sbnms.tiff",dpi=300,width=6,height=6,units="in")


###################################################################################

#total time period combined 2013-2018

#first get whale coordinates
sightings_xy <- data.frame(sightings_7to10$lon, sightings_7to10$lat)
coordinates(sightings_xy) <- c("sightings_7to10.lon", "sightings_7to10.lat")
proj4string(sightings_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
sightings_utm <- spTransform(sightings_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(sightings_utm, "SpatialPoints")

#bird stuff together
kps <- rbind(bw13orig, bw15orig, bw18orig)
kps <- filter(kps, spp=='GRSH')
kps <- kps %>% select(-spp)

kps$datetime <- strftime(strptime(kps$datetime,"%m/%d/%y %H:%M"),"%Y-%m-%d %H:%M:%S")

#format datetime as GMT
attr(kps$datetime, "tzone") <- 'GMT'

#make a local datetime (raw data is UTC)
kps$local_datetime <- kps$datetime
attr(kps$local_datetime, "tzone") <- 'EST'


all <- rbind(all_tracks_gom, kps)

#get coordinates for birds
all_xy <- data.frame(all$lon, all$lat)
coordinates(all_xy) <- c("all.lon", "all.lat")
proj4string(all_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
all_utm <- spTransform(all_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#estimate 90 and 50% UD for july2014 GOM data
all_kd <- kernelUD(all_utm,h="href", grid = 80, kern = c("bivnorm"), 
                          same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                          boundary = NULL)

#get the area of the UD
kernel.area(all_kd, percent = 50,
            unin = c("m"),
            unout = c("km2")) 

all_ud50 <- getverticeshr(all_kd,percent=50)



#figure out how many sightings fall within the all_tracks ud50
sightings_in_poly <- over(sightings_utm, all_ud50)

#find points outside polygon
points_outside_all_tracks<- sum(!complete.cases(sightings_in_poly))

#confirm that points outside is actually the real number
#nrow(sightings14_7to11_in_hulls[is.na(sightings14_7to11_in_hulls$id),])

#find percentage of sightings inside the all_tracks ud50
percent_overlap <- (nrow(sightings_7to10) - points_outside_all_tracks) / nrow(sightings_7to10) * 100

#plotting stuff

#convert to wgs84 for plotting
all_ud50_xy <- spTransform(all_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
all_ud50_fort <-fortify(all_ud50_xy)

##plot
t <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #whale sightings
  geom_point(data=sightings_7to10, aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) + 
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_ud50_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #scalebar
  ggsn::scalebar(data=all_ud50_fort,location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.05,st.size=3, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = all_ud50_fort, location = "bottomright", scale=0.15, symbol=1, anchor=c(x=-65.6, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="2013-2018") +
  annotate("text",x=-71.8,y=44.7,size=5, label = "95%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=12),
                                 plot.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

ggsave("./figures/total_sbnms.tiff",dpi=300,width=5,height=4,units="in")

##########################################################################################################

###weekly overlap

#2013

#huwh sightings 2013
sightings13 <- sightings_7to10 %>%
  filter(year=='2013') %>% droplevels()


#create this variable manually so code stays the same regardless of method for defining week
sightings13wk <- sightings13


#if using calendar weeks
sightings13wk<-sightings13wk %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#format week as date
sightings13wk$week <- as.Date(sightings13wk$week)

#Get week number of calendar year based on ISO8601
sightings13wk$week_id <- strftime(sightings13wk$week, format = "%V")

#bird tracks 2013 - refilter this to get the local time date(I put this in alltracks at beginning)
tracks13 <- filter(all, year=='2013') %>% droplevels()

#break up bird data by week

#create just a date column for tracks
tracks13$date <- as.Date(tracks13$local_datetime, format="%Y-%m-%d",tz="EST")

#for calendar week
tracks13wk<-tracks13 %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
tracks13wk$week_id <- strftime(tracks13wk$week, format = "%V")

#save full versions of both
sightings13wk_all <- sightings13wk
tracks13wk_all <- tracks13wk

#for overlap calcs, need to make sure same number of weeks in bird and sightings dataframe
#check weeks in both dataframes
unique(sightings13wk$week_id)
unique(tracks13wk$week_id)


sightings13wk <- sightings13wk %>%
  filter(week_id > 26) %>% droplevels()

tracks13wk <- tracks13wk %>%
  filter(week_id !=37)

#get the number of trackpoints in each week
tracks13wk %>% group_by(week_id) %>% summarise(trackptsperweek = n())

#get the number of huwh sightings in each week
sightings13wk %>% group_by(week_id) %>% summarise(sightperweek = n())

#let's chuck anything with less than 100 bird track points
tracks13wk <- tracks13wk %>%
  filter(week_id !=43)
sightings13wk <- sightings13wk %>%
  filter(week_id !=43)

#confirm matching weeks
unique(sightings13wk$week_id)
unique(tracks13wk$week_id)

###014
#2014

#huwh sightings 2014
sightings14 <- sightings_7to10 %>%
  filter(year=='2014') %>% droplevels()

#bird tracks 2014
tracks14 <- filter(all, year=='2014') %>% droplevels()

#create this variable manually so code stays the same regardless of method for defining week
sightings14wk <- sightings14

#if using calendar weeks
sightings14wk<-sightings14wk %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#format week as date
sightings14wk$week <- as.Date(sightings14wk$week)

#Get week number of calendar year based on ISO8601
sightings14wk$week_id <- strftime(sightings14wk$week, format = "%V")


#create just a date column for tracks
tracks14$date <- as.Date(tracks14$local_datetime, format="%Y-%m-%d",tz="EST")

#break by 7 days
#tracks14wk<-tracks14 %>% 
#mutate(week = as.factor(cut.Date(date, breaks="7 days")))

#for calendar week
tracks14wk<-tracks14 %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
tracks14wk$week_id <- strftime(tracks14wk$week, format = "%V")

#for overlap calcs, need to make sure same number of weeks in bird and sightings dataframe
#check weeks in both dataframes
unique(sightings14wk$week_id)
unique(tracks14wk$week_id)

#save full versions of both
sightings14wk_all <- sightings14wk
tracks14wk_all <- tracks14wk

#whales start on week 26 - week 43
#birds start week 29 - week 43
#put all non matching weeks as 0 in excel sheet
sightings14wk <- sightings14wk %>%
  filter(week_id > 28) %>% droplevels()

#confirm matching weeks
unique(sightings14wk$week_id)
unique(tracks14wk$week_id)

#get the number of trackpoints in each week
tracks14wk %>% group_by(week_id) %>% summarise(trackptsperweek = n())

#get the number of huwh sightings in each week
sightings14wk %>% group_by(week_id) %>% summarise(sightperweek = n())

#unique animals by week
tracks14wk %>% group_by(week_id) %>% summarise(count=n_distinct(tag))

#let's chuck anything less than 100 bird track points - no need for 2014



###2015

sightings15 <- sightings_7to10 %>%
  filter(year=='2015') %>% droplevels()

tracks15 <- filter(all, year=='2015') %>% droplevels()

sightings15wk <- sightings15

#if using calendar weeks
sightings15wk<-sightings15wk %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#format week as date
sightings15wk$week <- as.Date(sightings15wk$week)

#Get week number of calendar year based on ISO8601
sightings15wk$week_id <- strftime(sightings15wk$week, format = "%V")

#get the number of huwh sightings in each week
sightings15wk %>% group_by(week_id) %>% summarise(sightperweek = n())


#create just a date column for tracks
tracks15$date <- as.Date(tracks15$local_datetime, format="%Y-%m-%d",tz="EST")

#for calendar week
tracks15wk<-tracks15 %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
tracks15wk$week_id <- strftime(tracks15wk$week, format = "%V")

#for overlap calcs, need to make sure same number of weeks in bird and sightings dataframe
#check weeks in both dataframes
unique(sightings15wk$week_id)
unique(tracks15wk$week_id)

#save full versions of both
sightings15wk_all <- sightings15wk
tracks15wk_all <- tracks15wk
#whales start on week 26 - week 43
#birds start week 27 - week 43
#put all non matching weeks as 0 in excel sheet
sightings15wk <- sightings15wk %>%
  filter(week_id > 26) %>% droplevels()


#confirm matching weeks
unique(sightings15wk$week_id)
unique(tracks15wk$week_id)

#get the number of trackpoints in each week
tracks15wk %>% group_by(week_id) %>% summarise(trackptsperweek = n())

#unique animals by week
tracks15wk %>% group_by(week_id) %>% summarise(count=n_distinct(tag))

#chuck anything not over 100 pts - no need 2015


#2016

sightings16 <- sightings_7to10 %>%
  filter(year=='2016') %>% droplevels()

tracks16 <- filter(all, year=='2016') %>% droplevels()

sightings16wk <- sightings16

#if using calendar weeks
sightings16wk<-sightings16wk %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#format week as date
sightings16wk$week <- as.Date(sightings16wk$week)

#Get week number of calendar year based on ISO8601
sightings16wk$week_id <- strftime(sightings16wk$week, format = "%V")

#get the number of huwh sightings in each week
sightings16wk %>% group_by(week_id) %>% summarise(sightperweek = n())

#create just a date column for tracks
tracks16$date <- as.Date(tracks16$datetime, format="%Y-%m-%d",tz="EST")

#create just a date column for tracks
tracks16$date <- as.Date(tracks16$local_datetime, format="%Y-%m-%d",tz="EST")

#for calendar week
tracks16wk<-tracks16 %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
tracks16wk$week_id <- strftime(tracks16wk$week, format = "%V")

#for overlap calcs, need to make sure same number of weeks in bird and sightings dataframe
#check weeks in both dataframes
unique(sightings16wk$week_id)
unique(tracks16wk$week_id)

#save full versions of both
sightings16wk_all <- sightings16wk
tracks16wk_all <- tracks16wk

sightings16wk <- sightings16wk %>%
  filter(week_id > 26) %>% droplevels()
tracks16wk <- tracks16wk %>%
  filter(week_id < 42) %>% droplevels()

#confirm matching weeks
unique(sightings16wk$week_id)
unique(tracks16wk$week_id)

#get the number of trackpoints in each week
tracks16wk %>% group_by(week_id) %>% summarise(trackptsperweek = n())


###2017

##need to run this again because I remove July for the monthly stuff, but I want to count sightings by week just to have numbers
#filter sightings to july - nov to match bird tag time period
sightings_7to10 <- sightings %>%
  filter(month=="07" | month=="08" | month=="09" | month=="10")

sightings17 <- sightings_7to10 %>%
  filter(year=='2017') %>% droplevels()

tracks17 <- filter(all, year=='2017') %>% droplevels()

sightings17wk <- sightings17

#if using calendar weeks
sightings17wk<-sightings17wk %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#format week as date
sightings17wk$week <- as.Date(sightings17wk$week)

#Get week number of calendar year based on ISO8601
sightings17wk$week_id <- strftime(sightings17wk$week, format = "%V")

#get the number of huwh sightings in each week
sightings17wk %>% group_by(week_id) %>% summarise(sightperweek = n())

#create just a date column for tracks
tracks17$date <- as.Date(tracks17$local_datetime, format="%Y-%m-%d")

#create just a date column for tracks
tracks17$date <- as.Date(tracks17$local_datetime, format="%Y-%m-%d",tz="EST")

#for calendar week
tracks17wk<-tracks17 %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
tracks17wk$week_id <- strftime(tracks17wk$week, format = "%V")

#for overlap calcs, need to make sure same number of weeks in bird and sightings dataframe
#check weeks in both dataframes
unique(sightings17wk$week_id)
unique(tracks17wk$week_id)

#save full versions of both
sightings17wk_all <- sightings17wk
tracks17wk_all <- tracks17wk

sightings17wk <- sightings17wk %>%
  filter(week_id > 30) %>% droplevels()

#confirm matching weeks
unique(sightings17wk$week_id)
unique(tracks17wk$week_id)

#get the number of trackpoints in each week
tracks17wk %>% group_by(week_id) %>% summarise(trackptsperweek = n())



#2018

#huwh sightings 2013
sightings18 <- sightings_7to10 %>%
  filter(year=='2018') %>% droplevels()

#create this variable manually so code stays the same regardless of method for defining week
sightings18wk <- sightings18

#if using calendar weeks
sightings18wk<-sightings18wk %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#format week as date
sightings18wk$week <- as.Date(sightings18wk$week)

#Get week number of calendar year based on ISO8601
sightings18wk$week_id <- strftime(sightings18wk$week, format = "%V")

#bird tracks 2013 - refilter this to get the local time date(I put this in alltracks at beginning)
tracks18 <- filter(all, year=='2018') %>% droplevels()


#create just a date column for tracks
tracks18$date <- as.Date(tracks18$local_datetime, format="%Y-%m-%d",tz="EST")

#for calendar week
tracks18wk<-tracks18 %>% 
  mutate(week = as.factor(cut.Date(date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
tracks18wk$week_id <- strftime(tracks18wk$week, format = "%V")

#for overlap calcs, need to make sure same number of weeks in bird and sightings dataframe
#check weeks in both dataframes
unique(sightings18wk$week_id)
unique(tracks18wk$week_id)

#save full versions of both
sightings18wk_all <- sightings18wk
tracks18wk_all <- tracks18wk
#whales start on week 26 - week 43
#birds start week 29 - week 43
#put all non matching weeks as 0 in excel sheet
sightings18wk <- sightings18wk %>%
  filter(week_id > 26) %>% droplevels()

#confirm matching weeks
unique(sightings18wk$week_id)
unique(tracks18wk$week_id)

#get the number of trackpoints in each week
tracks18wk %>% group_by(week_id) %>% summarise(trackptsperweek = n())

#get the number of huwh sightings in each week
sightings18wk %>% group_by(week_id) %>% summarise(sightperweek = n())

#chuck anything with less than 100 bird points
tracks18wk <- tracks18wk %>%
  filter(week_id < 41) %>% droplevels()

sightings18wk <- sightings18wk %>%
  filter(week_id < 41) %>% droplevels()

#confirm matching weeks
unique(sightings18wk$week_id)
unique(tracks18wk$week_id)

##put all week dataframes together
sightings_wk <- as_tibble(rbind(sightings13wk,sightings14wk, sightings15wk, sightings16wk, sightings17wk,sightings18wk))
tracks_wk <- as_tibble(rbind(tracks13wk,tracks14wk, tracks15wk, tracks16wk, tracks17wk, tracks18wk))



#split dataframe by year and by week
dat_list_week <- group_split(sightings_wk, year, week)

list_xy_week <- lapply(dat_list_week, function(x){
  sightings_xy_week <- data.frame(x$lon, x$lat)
  coordinates(sightings_xy_week) <- c("x.lon", "x.lat")
  proj4string(sightings_xy_week) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  sightings_utm_week <- spTransform(sightings_xy_week, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(sightings_utm_week, "SpatialPoints")
})

bird_list_week <- group_split(tracks_wk, year, week)

bird_list_xy_week<- lapply(bird_list_week, function(y){
  bird_xy_week <- data.frame(y$lon, y$lat)
  coordinates(bird_xy_week) <- c("y.lon", "y.lat")
  proj4string(bird_xy_week) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  bird_utm_week <- spTransform(bird_xy_week, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(bird_utm_week, "SpatialPoints")
})

#take the new bird list and do all the kernel density stuff

all_tracks_kd_week <- lapply(bird_list_xy_week, function(y){
  kernelUD(y ,h="href", grid = 80, kern = c("bivnorm"), 
           same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
           boundary = NULL)#estimate 90 and 50% UD for july2014 GOM data
})

#find the area of each 50%UD
areas_week <- lapply(all_tracks_kd_week, function(y){ #figure out how to get month labels on here
  kernel.area(y, percent = 50,
              unin = c("m"),
              unout = c("km2"))
})

write.csv(areas_week, "areas_week_new.csv")

#take the new bird list and do all the kernel density stuff

all_tracks_kd_week <- lapply(bird_list_xy_week, function(y){
  kud_week <- kernelUD(y ,h="href", grid = 80, kern = c("bivnorm"), 
                       same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                       boundary = NULL) #estimate 90 and 50% UD for july2014 GOM data
  getverticeshr(kud_week,percent=50)#save the ud as a spatialpolygonsdataframe
})

#now I have a bird list with the hulls and a whale list in utm
#find sightings in the hulls

sightings_outside_week <- mapply(function(x,y){
  sightings_in_poly_week <- over(x,y)
  sightings_outside_poly_week <- sum(!complete.cases(sightings_in_poly_week))
},
list_xy_week, all_tracks_kd_week) 

#make overlap a dataframe -this is number of sightings outside the poly
sightings_outside_week <- data.frame(sightings_outside_week)
colnames(sightings_outside_week) <- c("points_outside")
sightings_outside_week <- split(sightings_outside_week,seq(nrow(sightings_outside_week))) #had to use seq(nrow) bc for some reason just split was dropping rows with 0 values

percent_overlap_week <- mapply(function(x,y){
  percent_overlap_week <- (nrow(x) - y$points_outside) / nrow(x) * 100
},
dat_list_week, sightings_outside_week)  

percent_overlap_week #get percent overlap for each year - percentage of  huwh sightings inside grsh 50ud

#save as a csv so I can copy and paste numbers correctly
write.csv(as.data.frame(percent_overlap_week), "percent_overlap_weeks_sbnms.csv")
write.csv(as.data.frame(sightings_outside_week), "points_outside_weeks_sbnms.csv")
###PLOTTING###

#convert to wgs84 for plotting
all_tracks_kd_xy_week <- lapply(all_tracks_kd_week, function(x){
  spTransform(x, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
})

#fortfy objects for ggplot
all_tracks_kd_xy_week_fort <- lapply(all_tracks_kd_xy_week, function(x){
  fortify(x)
})


#plots
#weekly
#2013
w1 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[1]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[1]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 27 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w2 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[2]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[2]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 28 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w3 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[3]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[3]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 29 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w4 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[4]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[4]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 30 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w5 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[5]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[5]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 31 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w6 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[6]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[6]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 32 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w7 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[7]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[7]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 33 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w8 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[8]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[8]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 34 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w9 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[9]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[9]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 35 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w10 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[10]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[10]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 36 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

#no whale sightings
w11 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[11]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[11]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 38 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w12 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[12]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[12]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 39 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[13]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[13]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 40 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[14]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[14]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 41 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[15]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[15]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 42 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



weeks13 <- ggarrange(w1 + rremove("ylab") + rremove("xlab"),
                     w2 + rremove("ylab") + rremove("xlab"),
                     w3 + rremove("ylab") + rremove("xlab"),
                     w4 + rremove("ylab") + rremove("xlab"),
                     w5 + rremove("ylab") + rremove("xlab"),
                     w6 + rremove("ylab") + rremove("xlab"),
                     w7 + rremove("ylab") + rremove("xlab"),
                     w8 + rremove("ylab") + rremove("xlab"),
                     w9 + rremove("ylab") + rremove("xlab"),
                     w10 + rremove("ylab") + rremove("xlab"),
                     w11 + rremove("ylab") + rremove("xlab"),
                     w12 + rremove("ylab") + rremove("xlab"),
                     w13 + rremove("ylab") + rremove("xlab"),
                     w14+ rremove("ylab") + rremove("xlab"),
                     w15 + rremove("ylab") + rremove("xlab"),
                    
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H","I","J","K","L","M","N","O"),
                     ncol=4,nrow=4) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(weeks13, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))

ggsave("./figures/2013_weekly_sbnms.tiff",dpi=300,width=10,height=8,units="in")



###2014 weeks

w1 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[16]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[16]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 29 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w2 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[17]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[17]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 30 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w3 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[18]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[18]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 31 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "98%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w4 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[19]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[19]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 32 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w5 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[20]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[20]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 33 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w6 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[21]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[21]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 34 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w7 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[22]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[22]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 35 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "97%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w8 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[23]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[23]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 36 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w9 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[24]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[24]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 37 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w10 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[25]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[25]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 38 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w11 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[26]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[26]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 39 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "63%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w12 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[27]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[27]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 40 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[28]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[28]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 41 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[29]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[29]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 42 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[30]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[30]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 43 2014") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

#multipanel

weeks14 <- ggarrange(w1 + rremove("ylab") + rremove("xlab"),
                     w2 + rremove("ylab") + rremove("xlab"),
                     w3 + rremove("ylab") + rremove("xlab"),
                     w4 + rremove("ylab") + rremove("xlab"),
                     w5 + rremove("ylab") + rremove("xlab"),
                     w6 + rremove("ylab") + rremove("xlab"),
                     w7 + rremove("ylab") + rremove("xlab"),
                     w8 + rremove("ylab") + rremove("xlab"),
                     w9 + rremove("ylab") + rremove("xlab"),
                     w10 + rremove("ylab") + rremove("xlab"),
                     w11 + rremove("ylab") + rremove("xlab"),
                     w12 + rremove("ylab") + rremove("xlab"),
                     w13 + rremove("ylab") + rremove("xlab"),
                     w14+ rremove("ylab") + rremove("xlab"),
                     w15 + rremove("ylab") + rremove("xlab"),
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H","I","J","K","L","M","N","O"),
                     ncol=4,nrow=4) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(weeks14, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))

ggsave("./figures/2014_weekly_sbnms.tiff",dpi=300,width=10,height=8,units="in")


###2015 weeks

w1 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[31]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[31]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 27 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))




w2 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[32]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[32]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 28 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "98%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w3 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[33]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[33]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 29 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w4 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[34]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[34]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 30 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w5 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[35]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[35]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 31 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "4%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w6 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[36]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[36]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 32 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "91%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w7 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[37]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[37]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 33 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "3%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w8 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[38]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[38]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 34 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w9 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[39]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[39]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 35 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w10 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[40]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[40]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 36 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "85%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w11 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[41]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[41]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 37 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w12 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[42]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[42]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 38 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[43]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[43]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 39 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[44]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[44]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 40 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[45]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[45]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 41 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "64%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[46]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[46]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 42 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w17 <- ggplot() +
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[47]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[47]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 43 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


#multipanel
weeks15 <- ggarrange(w1 + rremove("ylab") + rremove("xlab"),
                     w2 + rremove("ylab") + rremove("xlab"),
                     w3 + rremove("ylab") + rremove("xlab"),
                     w4 + rremove("ylab") + rremove("xlab"),
                     w5 + rremove("ylab") + rremove("xlab"),
                     w6 + rremove("ylab") + rremove("xlab"),
                     w7 + rremove("ylab") + rremove("xlab"),
                     w8 + rremove("ylab") + rremove("xlab"),
                     w9 + rremove("ylab") + rremove("xlab"),
                     w10 + rremove("ylab") + rremove("xlab"),
                     w11 + rremove("ylab") + rremove("xlab"),
                     w12 + rremove("ylab") + rremove("xlab"),
                     w13 + rremove("ylab") + rremove("xlab"),
                     w14+ rremove("ylab") + rremove("xlab"),
                     w15 + rremove("ylab") + rremove("xlab"),
                     w16 + rremove("ylab") + rremove("xlab"),
                     w17 + rremove("ylab") + rremove("xlab"),
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H","I","J","K","L","M","N","O","P","Q"),
                     ncol=4,nrow=5) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(weeks15, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))

ggsave("./figures/2015_weekly_sbnms.tiff",dpi=300,width=10,height=8,units="in")


###2016 weeks


w1 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[48]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[48]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 27 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "82%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w2 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[49]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[49]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 28 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w3 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[50]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[50]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 29 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w4 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[51]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[51]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 30 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "97%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w5 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[52]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[52]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 31 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "89%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w6 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[53]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[53]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 32 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "43%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w7 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #   #add in sightings
  geom_point(data=dat_list_week[[54]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[54]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 33 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w8 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[55]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[55]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 34 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "51%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w9 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[56]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[56]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 35 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w10 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[57]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[57]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 36 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w11 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[58]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[58]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 37 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w12 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #   #add in sightings
  geom_point(data=dat_list_week[[59]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[59]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 38 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[60]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[60]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 39 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "40%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[61]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[61]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 40 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[62]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[62]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="Week 41 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


#multipanel
weeks16 <- ggarrange(w1 + rremove("ylab") + rremove("xlab"),
                     w2 + rremove("ylab") + rremove("xlab"),
                     w3 + rremove("ylab") + rremove("xlab"),
                     w4 + rremove("ylab") + rremove("xlab"),
                     w5 + rremove("ylab") + rremove("xlab"),
                     w6 + rremove("ylab") + rremove("xlab"),
                     w7 + rremove("ylab") + rremove("xlab"),
                     w8 + rremove("ylab") + rremove("xlab"),
                     w9 + rremove("ylab") + rremove("xlab"),
                     w10 + rremove("ylab") + rremove("xlab"),
                     w11 + rremove("ylab") + rremove("xlab"),
                     w12 + rremove("ylab") + rremove("xlab"),
                     w13 + rremove("ylab") + rremove("xlab"),
                     w14+ rremove("ylab") + rremove("xlab"),
                     w15 + rremove("ylab") + rremove("xlab"),
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H","I","J","K","L","M","N","O"),
                     ncol=4,nrow=4) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(weeks16, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))

ggsave("./figures/2016_weekly_sbnms.tiff",dpi=300,width=10,height=8,units="in")



###2017 weeks



w1 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[63]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[63]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 31 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "50%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w2 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[64]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[64]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 32 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w3 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[65]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[65]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 33 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w4 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[66]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[66]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 34 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "95%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w5 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[67]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[67]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 35 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "64%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w6 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[68]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[68]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 36 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "12%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w7 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[69]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[69]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 37 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w8 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[70]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[70]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 38 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "92%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w9 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[71]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[71]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 39 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "69%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w10 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[72]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[72]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
 
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 40 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "68%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w11 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[73]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[73]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 41 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w12 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[74]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[74]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 42 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[75]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[75]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="Week 43 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

#multipanel
weeks17 <- ggarrange(w1 + rremove("ylab") + rremove("xlab"),
                     w2 + rremove("ylab") + rremove("xlab"),
                     w3 + rremove("ylab") + rremove("xlab"),
                     w4 + rremove("ylab") + rremove("xlab"),
                     w5 + rremove("ylab") + rremove("xlab"),
                     w6 + rremove("ylab") + rremove("xlab"),
                     w7 + rremove("ylab") + rremove("xlab"),
                     w8 + rremove("ylab") + rremove("xlab"),
                     w9 + rremove("ylab") + rremove("xlab"),
                     w10 + rremove("ylab") + rremove("xlab"),
                     w11 + rremove("ylab") + rremove("xlab"),
                     w12 + rremove("ylab") + rremove("xlab"),
                     w13 + rremove("ylab") + rremove("xlab"),
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H","I","J","K","L","M"),
                     ncol=4,nrow=4) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(weeks17, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))

ggsave("./figures/2017_weekly_sbnms.tiff",dpi=300,width=10,height=8,units="in")


###2018 weeks
world_map <- map_data("world")

w1 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[76]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[76]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 27 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))



w2 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[77]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[77]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 28 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w3 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[78]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[78]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 29 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w4 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[79]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[79]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 30 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "99%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w5 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[80]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[80]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 31 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w6 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[81]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[81]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="Week 32 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w7 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[82]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[82]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="Week 33 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w8 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[83]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[83]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 34 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w9 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[84]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[84]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 35 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w10 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[85]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[85]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +

  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 36 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

w11 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[86]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[86]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="Week 37 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "7%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w12 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[87]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[87]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="Week 38 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "0%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[88]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[88]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="Week 39 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "97%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


w14 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=dat_list_week[[89]], aes(x=lon, y=lat),inherit.aes = FALSE,
             size=0.1, show.legend=FALSE) +
  #ud50 for all tracks combined in GOM
  geom_polygon(data = all_tracks_kd_xy_week_fort[[89]],
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='dodgerblue2', alpha=0.25) + 
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="Week 40 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "71%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))





#multipanel
weeks18 <- ggarrange(w1 + rremove("ylab") + rremove("xlab"),
                     w2 + rremove("ylab") + rremove("xlab"),
                     w3 + rremove("ylab") + rremove("xlab"),
                     w4 + rremove("ylab") + rremove("xlab"),
                     w5 + rremove("ylab") + rremove("xlab"),
                     w6 + rremove("ylab") + rremove("xlab"),
                     w7 + rremove("ylab") + rremove("xlab"),
                     w8 + rremove("ylab") + rremove("xlab"),
                     w9 + rremove("ylab") + rremove("xlab"),
                     w10 + rremove("ylab") + rremove("xlab"),
                     w11 + rremove("ylab") + rremove("xlab"),
                     w12 + rremove("ylab") + rremove("xlab"),
                     w13 + rremove("ylab") + rremove("xlab"),
                     w14+ rremove("ylab") + rremove("xlab"),
                   
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H","I","J","K","L","M","N"),
                     ncol=4,nrow=4) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(weeks18, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))

ggsave("./figures/2018_weekly_sbnms.tiff",dpi=300,width=10,height=8,units="in")




