#This script uses great shearwater satellite telemetry data from SBNMS and humpback whale sightings from CCS
#to: calculate 50% utilization distributions for shearwaters within the Gulf of Maine from July - October, 
#calculate how many and the percentage of humpback sightings from CCS GOM surveysthat fall within the 50UD
#Tammy Silva 10/12/2021


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

#get bathy from marmap so I can have contours
bath_cont <- getNOAA.bathy (lon1 = -65, lon2 = -72.5, lat1 = 40, lat2 = 45, resolution = 1)
bath_cont_fort<- fortify.bathy(bath_cont) #for ggplot

#read in 40m contou shapefile
#bath40<-readOGR("./data/bathy40m.shp") 
#bath40 <- spTransform(bath40, proj_wgs84)
#bath40_fort <- fortify(bath40) #for ggplot

##################################################################################

#load in GOM survey effort data polygons made by Mike T
eff13<-readOGR("./data/buffer_10km_2013.shp") 
eff13 <- spTransform(eff13, proj_wgs84)
eff13_fort <- fortify(eff13) #for ggplot

eff14<-readOGR("./data/buffer_10km_2014.shp") 
eff14 <- spTransform(eff14, proj_wgs84)
eff14_fort <- fortify(eff14) #for ggplot

eff15<-readOGR("./data/buffer_10km_2015.shp") 
eff15 <- spTransform(eff15, proj_wgs84)
eff15_fort <- fortify(eff15) #for ggplot

eff16<-readOGR("./data/buffer_10km_2016.shp") 
eff16 <- spTransform(eff16, proj_wgs84)
eff16_fort <- fortify(eff16) #for ggplot

eff17<-readOGR("./data/buffer_10km_2017.shp") 
eff17 <- spTransform(eff17, proj_wgs84)
eff17_fort <- fortify(eff17) #for ggplot

eff18<-readOGR("./data/buffer_10km_2018.shp") 
eff18 <- spTransform(eff18, proj_wgs84)
eff18_fort <- fortify(eff18) #for ggplot

#####################################################################
###whale stuff

###whale stuff
### load in CCS sightings data
sightings <- read_excel("./data/CCS GOM SURVEY DATA FOR DOM.xlsx",na=" ") 

#formatting
sightings <- rename(sightings, year = Year, month = Month, day = Day, lat = Latitude, lon = Longitude )


#make a date column
sightings <- mutate(sightings,date = paste(year, month, day, sep = "-"))

#format date as dates
sightings$date <- as.Date(sightings$date,"%Y-%m-%d")

#check summary
summary(sightings)

#drop one row with NA for position
sightings <- drop_na(sightings)

#format month as month
sightings$month <- as.factor(sightings$month)

#filter sightings to july - nov to match bird tag time period
sightings_7to10 <- sightings %>%
  filter(month=="7" | month=="8" | month=="9" | month=="10")

#get some sighting sumary information
sightings_7to10 %>% group_by(year) %>% summarise(sightperyear = n())
numsight<- sightings_7to10 %>% group_by(year, month) %>% summarise(sightperyearmonth = n())



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



#bird stuff

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
#bw13 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())

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
sight13jul<- subset(sight13, sight13$month == '7')

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
points13jul_outside_all_tracks <- sum(!complete.cases(sight13jul_in_poly))

points13jul_outside_all_tracks

#find percentage of sightings inside the all_tracks ud50
percent_overlap13jul <- (nrow(sight13jul) - points13jul_outside_all_tracks) / nrow(sight13jul_in_poly) * 100

percent_overlap13jul

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
sight13aug<- subset(sight13, sight13$month == '8')

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

points13aug_outside_all_tracks

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
sight13sep<- subset(sight13, sight13$month == '9')

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

points13sep_outside_all_tracks

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

points13oct_outside_all_tracks

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

points15_outside_all_tracks

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


###months
#2015

#July 2015

#whale data

#subset whale sightings
sight15jul<- subset(sight15, sight15$month == '7')

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
points15jul_outside_all_tracks

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
sight15aug<- subset(sight15, sight15$month == '8')

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
sight15sep<- subset(sight15, sight15$month == '9')

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
sight18jul<- subset(sight18, sight18$month == '7')

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
sight18aug<- subset(sight18, sight18$month == '8')

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
sight18sep<- subset(sight18, sight18$month == '9')

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


##2014-2017
#years and months

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
  filter(!(year=='2017' & month=='7')) %>%
  droplevels()

#no sightings in ocober 2014 -  remove birds from this month
all_tracks_gom <- all_tracks_gom %>%
  filter(!(year=='2014' & month=='10')) %>%
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
#write.csv(areas, "areas_months14to17GOM.csv") 

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
  
  #effort
  geom_polygon(data = eff13_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='gray47', alpha=0.25) + 
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
  annotate("text",x=-71,y=44.6,size=5, label = "80%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  


#2014

p14 <-  ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #effort
  geom_polygon(data = eff14_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='gray47', alpha=0.25) + 
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
  annotate("text",x=-71,y=44.6,size=5, label = "80%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 



#2015

p15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #effort
  geom_polygon(data = eff15_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='gray47', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "88%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 


#2016

p16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #effort
  geom_polygon(data = eff16_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='gray47', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "87%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 


# 2017

p17 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #effort
  geom_polygon(data = eff17_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='gray47', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "77%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))   


# 2018

p18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  
  #effort
  geom_polygon(data = eff18_fort,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='gray47', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "26%",color="white") +
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
ggsave("./figures/all_years13to18_gom.tiff",dpi=300,width=6,height=8.5,units="in")

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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +

  labs(x="Longitude", y="Latitude",title="July 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "71%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "34%",color="white") +
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
  #bathym40m
 # geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
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
ggsave("./figures/2013_monthly_gom.tiff",dpi=300,width=6,height=6,units="in")

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
  #bathym40m
  #geom_path(data = bath40_fort,
           # aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
 
  labs(x="Longitude", y="Latitude",title="July 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "79%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "94%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
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
ggsave("./figures/2015_monthly_gom.tiff",dpi=300,width=6,height=6,units="in")


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
  #bathym40m
  #geom_path(data = bath40_fort,
           # aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title=" July 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "55%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
           # aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "5%",color="white") +
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
  #bathym40m
  #geom_path(data = bath40_fort,
            #aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
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
  #bathym40m
 # geom_path(data = bath40_fort,
           # aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'black', size = .2) +
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
ggsave("./figures/2018_monthly_gom.tiff",dpi=300,width=6,height=6,units="in")

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
  annotate("text",x=-71,y=44.6,size=5, label = "53%",color="white") +
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
  annotate("text",x=-71,y=44.6,size=5, label = "9%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months14gom <- ggarrange(j14 + rremove("ylab") + rremove("xlab"),
                      a14+ rremove("ylab") + rremove("xlab"),
                      s14+ rremove("ylab") + rremove("xlab"),
                    
                      
                      labels = c("A", "B", "C"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months14gom, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2014_monthly_gom.tiff",dpi=300,width=6,height=6,units="in")



###2016

j16 <- ggplot() + 
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
  #scalebar
  ggsn::scalebar(data=all_tracks_kd_xy_month_fort[[4]],location="bottomright", st.bottom=TRUE,dist=50, dist_unit="km",height=0.03,
                 st.dist = 0.15,st.size=2, box.fill=c("black", "white"), transform=TRUE, model="WGS84", anchor=c(x=-65.4, y=40.5)) + 
  north(data = all_tracks_kd_xy_month_fort[[4]], location = "bottomright", scale=0.6, symbol=1, anchor=c(x=-65.5, y=40.7)) +
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="July 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

###Plot Aug 2016

a16 <- ggplot() + 
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
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="August 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "90%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Sept 2016

s16 <- ggplot() + 
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
  
  labs(x="Longitude", y="Latitude",title="September 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "53%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot oct 2016

o16 <- ggplot() + 
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
  
  labs(x="Longitude", y="Latitude",title="October 2016") +
  annotate("text",x=-71,y=44.6,size=5, label = "33%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months16gom <- ggarrange(j16 + rremove("ylab") + rremove("xlab"),
                      a16+ rremove("ylab") + rremove("xlab"),
                      s16+ rremove("ylab") + rremove("xlab"),
                      o16+ rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months16gom, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2016_monthly_gom.tiff",dpi=300,width=6,height=6,units="in")


###Plot for July 2017 - just sightings
#refilter to get sightings
#sight717 <- sightings %>% filter(year=='2017' & month=='07')

#j17 <- ggplot() + 
  
  #land
  # geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  # coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  # #bathym40m
  # geom_path(data = bath40_fort,
  #           aes(x = long, y = lat, group = group),inherit.aes = FALSE, color = 'darkgray', size = .2) +
  # #sbnms
  # geom_path(data = sbnms_fort,
  #           aes(x = long, y = lat, group = group),inherit.aes = FALSE,
  #           color = 'black', size = .2) +
  # #   #add in sightings
  # geom_point(data=sight717, aes(x=lon, y=lat),inherit.aes = FALSE,
  #            size=0.1, show.legend=FALSE) +
  # labs(x="Longitude", y="Latitude",title="July 2017") +
  # annotate("text",x=-70.3,y=44.6,size=5, label = "No tagged birds",color="white") +
  # theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
  #                                plot.title=element_text(size=12))


###Plot Aug 2017

a17 <- ggplot() + 
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
  annotate("text",x=-71,y=44.6,size=5, label = "51%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Sept 2017

s17 <- ggplot() + 
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
  #land
  geom_polygon(data=world_map,aes(x = long, y = lat, group=group), fill="gray47") + 
  coord_sf(xlim = c(-72.5, -65), ylim = c(40, 45), expand = FALSE) +
  #sbnms
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="September 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))


###Plot Oct 2017
o17 <- ggplot() + 
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
  
  labs(x="Longitude", y="Latitude",title="October 2017") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months17gom <- ggarrange(
                      a17+ rremove("ylab") + rremove("xlab"),
                      s17+ rremove("ylab") + rremove("xlab"),
                      o17+ rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months17gom, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2017_monthly_gom.tiff",dpi=300,width=6,height=6,units="in")


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

#need to run tracks code again bc I filtered Oct 2014 out of the data bc I had no whale sighings
#need to put it back in

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
points_outside_all_tracks
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
  
  geom_path(data = sbnms_fort,
            aes(x = long, y = lat, group = group),inherit.aes = FALSE,
            color = 'black', size = .2) +
  
  labs(x="Longitude", y="Latitude",title="2013-2018") +
  annotate("text",x=-71.8,y=44.7,size=5, label = "86%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=12),
                                 plot.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

ggsave("./figures/total_gom.tiff",dpi=300,width=5,height=4,units="in")


############################################################################################

#plot HUWH entanglements
#data from NEFSC / DWiley

ent <- read_excel("./data/2013-2018 HI GoM HUWH for D.Wiley.xlsx",na=" ") 

str(ent)

#get rid of notes cols for now
en <- ent %>%
  select(Date, ID, Location, Latitude, Longitude, Cause, Fate)

#change col names
en <- en %>%
  rename(lat = Latitude, lon = Longitude, date = Date)

#fill in Stellwagen Bank lat lon with generic middle of Bank coordinates
#using coordinates for sand lance site C3 42.278	70.28

en <- en %>%
  mutate(lat = ifelse(grepl("Stellwagen", Location),"42.278", lat))
en <- en %>%
  mutate(lon = ifelse(grepl("Stellwagen", Location),"70.28", lon))

#lat lon as numeric
en$lat <- as.numeric(en$lat)
en$lon <- as.numeric(en$lon)

str(en)

#create year and month week cols
en$Date <- as.Date(en$date) #this is to make the cut.Date work below
en$year <- format(en$Date, "%Y")
en$month <- format(en$Date, "%m")

#Month as factor
en$month <- as.factor(en$month)

#cause as factor
en$Cause <- as.factor(en$Cause)

#filter out just use entanglements
en <- en %>%
  filter(Cause=='EN')

#limit to GOM
en <- filter(en, lat > 40 & lat < 45) 
en <- filter(en, lon < -65)

#for calendar week
en <- en %>% 
  mutate(week = as.factor(cut.Date(Date, breaks="week", start.on.monday = FALSE)))

#Get week number of calendar year based on ISO8601
en$week_id <- strftime(en$week , format = "%V")

#get ent from july to october to match sightings
en <- en %>%
  filter(month=='07' |month=='08' | month=='09' | month=='10' )

#convert ent locations to utm to calculate overlap
#first get ent coordinates
en_xy <- data.frame(en$lon, en$lat)
coordinates(en_xy) <- c("en.lon", "en.lat")
proj4string(en_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
en_utm <- spTransform(en_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(en_utm, "SpatialPoints")

#get overlap between total grsh 50 ud and entanglements
en_in_poly <- over(en_utm, all_ud50)

#find points outside polygon
en_outside_all_tracks<- sum(!complete.cases(en_in_poly))

#confirm that points outside is actually the real number
#nrow(sightings14_7to11_in_hulls[is.na(sightings14_7to11_in_hulls$id),])

#find percentage of sightings inside the all_tracks ud50
percent_overlap_en <- (nrow(en) - en_outside_all_tracks) / nrow(en) * 100

#sbnms to utm
sbnms_utm <- spTransform(sbnms, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

#get number of ent inside SBNMS
en_in_sbnms <- over(en_utm, sbnms_utm)

#find points outside polygon
en_outside_sbnms<- sum(!complete.cases(en_in_sbnms))




##plot total ent
ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #add HI locations
  geom_point(data=en, aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  
  labs(x="Longitude", y="Latitude",title="Observed Entanglements 2013 - 2018") +
  annotate("text",x=-71,y=44.6,size=5, label = "63%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

ggsave("./figures/total_entanglements.tiff",width=4,height=4,units="in")

#en by year
en13 <- en %>% filter(year=='2013')
en14 <- en %>% filter(year=='2014')
en15 <- en %>% filter(year=='2015')
en16 <- en %>% filter(year=='2016')
en17 <- en %>% filter(year=='2017')
en18 <- en %>% filter(year=='2018')

#en %>% group_by(year) %>% summarise(entperyear = n())

#make points utm
#2013
en13_xy <- data.frame(en13$lon, en13$lat)
coordinates(en13_xy) <- c("en13.lon", "en13.lat")
proj4string(en13_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
en13_utm <- spTransform(en13_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(en13_utm, "SpatialPoints")

#get overlap
en13_in_poly <- over(en13_utm, bw13_ud50gs)

#find points outside polygon
en13_outside_all_tracks<- sum(!complete.cases(en13_in_poly))

#find percentage of sightings inside the all_tracks ud50
percent_overlapen13 <- (nrow(en13) - en13_outside_all_tracks) / nrow(en13_in_poly) * 100

#2015
en15_xy <- data.frame(en15$lon, en15$lat)
coordinates(en15_xy) <- c("en15.lon", "en15.lat")
proj4string(en15_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
en15_utm <- spTransform(en15_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(en15_utm, "SpatialPoints")

#get overlap
en15_in_poly <- over(en15_utm, bw15_ud50gs)

#find points outside polygon
en15_outside_all_tracks<- sum(!complete.cases(en15_in_poly))

#find percentage of sightings inside the all_tracks ud50
percent_overlapen15 <- (nrow(en15) - en15_outside_all_tracks) / nrow(en15_in_poly) * 100

#2018
en18_xy <- data.frame(en18$lon, en18$lat)
coordinates(en18_xy) <- c("en18.lon", "en18.lat")
proj4string(en18_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#transform to utm
en18_utm <- spTransform(en18_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))

#create spatial points object
as(en18_utm, "SpatialPoints")

#get overlap
en18_in_poly <- over(en18_utm, bw18_ud50gs)

#find points outside polygon
en18_outside_all_tracks<- sum(!complete.cases(en18_in_poly))

#find percentage of sightings inside the all_tracks ud50
percent_overlapen18 <- (nrow(en18) - en18_outside_all_tracks) / nrow(en18_in_poly) * 100


#2014-2017
en14to17 <- en %>% filter(year=='2014' | year=='2016' | year=='2017') %>% droplevels()
ent_list <- split(en14to17, en14to17$year) 
list_xy<- lapply(ent_list, function(x){
  sightings_xy <- data.frame(x$lon, x$lat)
  coordinates(sightings_xy) <- c("x.lon", "x.lat")
  proj4string(sightings_xy) <- CRS(" +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  sightings_utm <- spTransform(sightings_xy, CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  as(sightings_utm, "SpatialPoints")
})



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


#plot entanglements for each year
e13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=en13, aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  annotate("text",x=-71,y=44.6,size=5, label = "40%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  




world_map <- map_data("world")

e15 <-  ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=en15, aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  annotate("text",x=-71,y=44.6,size=5, label = "71%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 





e18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=en18, aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  annotate("text",x=-71,y=44.6,size=5, label = "13%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 



e14 <- ggplot() +
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=ent_list[[1]], aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  annotate("text",x=-71,y=44.6,size=5, label = "33%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 




e16 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=ent_list[[2]], aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  annotate("text",x=-71,y=44.6,size=5, label = "81%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))   




e17 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #   #add in sightings
  geom_point(data=ent_list[[3]], aes(x=lon, y=lat),color="red",inherit.aes = FALSE,
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
  annotate("text",x=-71,y=44.6,size=5, label = "70%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 

#make multipanel plot of annual overlap
#amazing reference / help for all you can do in ggarrange 
#https://www.r-bloggers.com/2017/07/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

ent_all_years <- ggarrange(e13 + rremove("ylab") + rremove("xlab"),
                           e14+ rremove("ylab") + rremove("xlab"),
                           e15+ rremove("ylab") + rremove("xlab"),
                           e16+ rremove("ylab") + rremove("xlab"),
                           e17+ rremove("ylab") + rremove("xlab"),
                           e18 + rremove("ylab") + rremove("xlab"),
                           labels = c("A", "B", "C", "D", "E", "F"),
                           ncol=2,nrow=3) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(ent_all_years, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/ent13to18.tiff",dpi=300,width=6,height=8.5,units="in")
