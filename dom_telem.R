#this file recreates the DOM analysis by KP and quantifies spatial overlap between GRSH and HUWH using satellite telemetry
#data from 2013, 2015, 2018
#Tammy Silva
#11/3/2021

library(adehabitatHR) #can export shape files of polygons
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(tidyverse)
library(readxl)
library(ggpubr)
library(marmap)

#projections
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


##some summary values

#get trackpts and idivdiuals by year
bw13 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animalsby year
bw13 %>% group_by(spp,year) %>% summarise(count=n_distinct(tag))
# 
# # #get track pts for each month
bw13 %>% group_by(spp,year, month) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animals by month
bw13 %>% group_by(spp,year, month) %>% summarise(count=n_distinct(tag))


#get locations ready for UD calc

coordinates(bw13)<- ~x+y		# set the projection on UTM coords
proj4string(bw13)<- proj_utm

#calculae kernel UD - using settings developed by KP
kudl13 <- kernelUD(bw13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                 same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                 boundary = NULL)


kernel.area(kudl13, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudl13, method = c("HR"), percent = 50, conditional = FALSE)

bw13_ud50 <- getverticeshr(kudl13, percent = 50)

#for plotting
#convert to wgs84 for plotting
bw13_ud50_xy <- spTransform(bw13_ud50, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# 
# #sightings_xy is already in wgs84
# 
# #fortify objects for ggplot
bw13_ud50_fort <-fortify(bw13_ud50_xy)

#need separate data bc of this stuuuuupid code
gs13 <-filter(bw13_ud50_fort, id=="GRSH")
hw13 <-filter(bw13_ud50_fort, id=="HUWH")


###July 2013
# July
jul13<- subset(bw13, bw13$month == 7)

#calculae kernel UD - using settings developed by KP
kudljul13 <- kernelUD(jul13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                   same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                   boundary = NULL)


kernel.area(kudljul13, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudljul13, method = c("HR"), percent = 50, conditional = FALSE)

jul13_ud50 <- getverticeshr(kudljul13, percent = 50)

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
hwjul13 <-filter(jul13_ud50_fort, id=="HUWH")



###aug 2013

aug13<- subset(bw13, bw13$month == 8)

#calculae kernel UD - using settings developed by KP
kudlaug13 <- kernelUD(aug13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudlaug13, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudlaug13, method = c("HR"), percent = 50, conditional = FALSE)

aug13_ud50 <- getverticeshr(kudlaug13, percent = 50)

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
hwaug13 <-filter(aug13_ud50_fort, id=="HUWH")


###sep 2013

sep13<- subset(bw13, bw13$month == 9)

#calculae kernel UD - using settings developed by KP
kudlsep13 <- kernelUD(sep13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudlsep13, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudlsep13, method = c("HR"), percent = 50, conditional = FALSE)

sep13_ud50 <- getverticeshr(kudlsep13, percent = 50)

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
hwsep13 <-filter(sep13_ud50_fort, id=="HUWH")


###oct 2013
#oct
oct13<- subset(bw13, bw13$month == 10)

#calculae kernel UD - using settings developed by KP
kudloct13 <- kernelUD(oct13[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudloct13, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudloct13, method = c("HR"), percent = 50, conditional = FALSE)

oct13_ud50 <- getverticeshr(kudloct13, percent = 50)

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
hwoct13 <-filter(oct13_ud50_fort, id=="HUWH")



###2015

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


##some summary values

#get trackpts and idivdiuals by year
 bw15 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animalsby year
bw15 %>% group_by(spp,year) %>% summarise(count=n_distinct(tag))
# 
# # #get track pts for each month
bw15 %>% group_by(spp,year, month) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animals by month
 bw15 %>% group_by(spp,year, month) %>% summarise(count=n_distinct(tag))


#get locations ready for UD calc

coordinates(bw15)<- ~x+y		# set the projection on UTM coords
proj4string(bw15)<- proj_utm

#calculae kernel UD - using settings developed by KP
kudl15 <- kernelUD(bw15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                   same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                   boundary = NULL)


kernel.area(kudl15, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudl15, method = c("HR"), percent = 50, conditional = FALSE)

bw15_ud50 <- getverticeshr(kudl15, percent = 50)

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


###July 2015

jul15<- subset(bw15, bw15$month == 7)

#calculae kernel UD - using settings developed by KP
kudljul15 <- kernelUD(jul15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudljul15, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudljul15, method = c("HR"), percent = 50, conditional = FALSE)

jul15_ud50 <- getverticeshr(kudljul15, percent = 50)

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
hwjul15 <-filter(jul15_ud50_fort, id=="HUWH")


###aug 2015

aug15<- subset(bw15, bw15$month == 8)

#calculae kernel UD - using settings developed by KP
kudlaug15 <- kernelUD(aug15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudlaug15, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudlaug15, method = c("HR"), percent = 50, conditional = FALSE)

aug15_ud50 <- getverticeshr(kudlaug15, percent = 50)

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
hwaug15 <-filter(aug15_ud50_fort, id=="HUWH")



###sep 2015

sep15<- subset(bw15, bw15$month == 9)

#calculae kernel UD - using settings developed by KP
kudlsep15 <- kernelUD(sep15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudlsep15, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudlsep15, method = c("HR"), percent = 50, conditional = FALSE)

sep15_ud50 <- getverticeshr(kudlsep15, percent = 50)

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
hwsep15 <-filter(sep15_ud50_fort, id=="HUWH")



###oct 2015

oct15<- subset(bw15, bw15$month == 10)

#calculae kernel UD - using settings developed by KP
kudloct15 <- kernelUD(oct15[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudloct15, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudloct15, method = c("HR"), percent = 50, conditional = FALSE)

oct15_ud50 <- getverticeshr(kudloct15, percent = 50)

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
hwoct15 <-filter(oct15_ud50_fort, id=="HUWH")


###2018

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

#need to remove 4 whales tagged in BOF that may show bias
 bw18 <- filter(bw18, tag !='EGREQUE' & tag !='SUTURES' & tag !='GODZILLA', tag !='PUPPET' )


##some summary values

#get trackpts and idivdiuals by year
bw18 %>% group_by(spp,year) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animalsby year
bw18 %>% group_by(spp,year) %>% summarise(count=n_distinct(tag))
# 
# # #get track pts for each month
bw18 %>% group_by(spp,year, month) %>% summarise(ptsperyearmonth = n())
# 
# 
# #unique animals by month
bw18 %>% group_by(spp,year, month) %>% summarise(count=n_distinct(tag))


#get locations ready for UD calc

coordinates(bw18)<- ~x+y		# set the projection on UTM coords
proj4string(bw18)<- proj_utm

#calculae kernel UD - using settings developed by KP
kudl18 <- kernelUD(bw18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                   same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                   boundary = NULL)


kernel.area(kudl18, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudl18, method = c("HR"), percent = 50, conditional = FALSE)

bw18_ud50 <- getverticeshr(kudl18, percent = 50)

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
hw18 <-filter(bw18_ud50_fort, id=="HUWH")


###July 2018

jul18<- subset(bw18, bw18$month == 7)

#calculae kernel UD - using settings developed by KP
kudljul18 <- kernelUD(jul18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudljul18, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudljul18, method = c("HR"), percent = 50, conditional = FALSE)

jul18_ud50 <- getverticeshr(kudljul18, percent = 50)

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
hwjul18 <-filter(jul18_ud50_fort, id=="HUWH")


###aug 2018

aug18<- subset(bw18, bw18$month == 8)

#calculae kernel UD - using settings developed by KP
kudlaug18 <- kernelUD(aug18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudlaug18, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudlaug18, method = c("HR"), percent = 50, conditional = FALSE)

aug18_ud50 <- getverticeshr(kudlaug18, percent = 50)

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
hwaug18 <-filter(aug18_ud50_fort, id=="HUWH")


###sep 2018

sep18<- subset(bw18, bw18$month == 9)

#calculae kernel UD - using settings developed by KP
kudlsep18 <- kernelUD(sep18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudlsep18, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudlsep18, method = c("HR"), percent = 50, conditional = FALSE)

sep18_ud50 <- getverticeshr(kudlsep18, percent = 50)

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
hwsep18 <-filter(sep18_ud50_fort, id=="HUWH")


###oct 2018

oct18<- subset(bw18, bw18$month == 10)

#calculae kernel UD - using settings developed by KP
kudloct18 <- kernelUD(oct18[,2], h="href", grid = 80, kern = c("bivnorm"), 
                      same4all = TRUE, hlim = c(0.1, 1.5), extent = 1, 
                      boundary = NULL)


kernel.area(kudloct18, percent = 50, 
            unin = c("m"),
            unout = c("km2"))

kerneloverlaphr(kudloct18, method = c("HR"), percent = 50, conditional = FALSE)

oct18_ud50 <- getverticeshr(kudloct18, percent = 50)

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
hwoct18 <-filter(oct18_ud50_fort, id=="HUWH")


###################################################################
###plot years

### 2013 

world_map <- map_data("world")

p13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hw13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 


###2015

p15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hw15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 

###2018 

p18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hw18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "15%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12)) 


years13telem <- ggarrange(p13 + rremove("ylab") + rremove("xlab"),
                      p15 + rremove("ylab") + rremove("xlab"),
                      p18 + rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(years13telem, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/all_years_telem.tiff",dpi=300,width=6,height=6,units="in")

##########################################################################

##monthly plots

##2013

###plot july 2013

world_map <- map_data("world")

jul13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwjul13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  labs(x="Longitude", y="Latitude",title=" July 2013") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  

### august 2013 



aug13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwaug13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  


###sept 2013 

world_map <- map_data("world")

sep13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwsep13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

###oct2013 

world_map <- map_data("world")

oct13 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwoct13,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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

months13telem <- ggarrange(jul13 + rremove("ylab") + rremove("xlab"),
                           aug13 + rremove("ylab") + rremove("xlab"),
                           sep13 + rremove("ylab") + rremove("xlab"),
                           oct13 + rremove("ylab") + rremove("xlab"),
                           
                           labels = c("A", "B", "C", "D"),
                           ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months13telem, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2013_monthly_overlap_telem.tiff",dpi=300,width=6,height=6,units="in")



###2015

###july 2015

jul15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwjul15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  labs(x="Longitude", y="Latitude",title=" July 2015") +
  annotate("text",x=-71,y=44.6,size=5, label = "100%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  

###august 2015 

aug15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwaug15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "88%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  


###sept 2015 

sep15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwsep15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "95%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

### oct2015 

world_map <- map_data("world")

oct15 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwoct15,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "73%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

months15telem <- ggarrange(jul15 + rremove("ylab") + rremove("xlab"),
                      aug15 + rremove("ylab") + rremove("xlab"),
                      sep15 + rremove("ylab") + rremove("xlab"),
                      oct15 + rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months15telem, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2015_monthly_overlap_telem.tiff",dpi=300,width=6,height=6,units="in")


###2018
###july 2018


jul18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwjul18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "12%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))

###august 2018 

aug18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwaug18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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
  annotate("text",x=-71,y=44.6,size=5, label = "10%",color="white") +
  theme_bw(base_size=12) + theme(axis.text=element_text(size=10),
                                 plot.title=element_text(size=12))  


###sept 2018 

sep18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwsep18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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

###oct2018 


oct18 <- ggplot() + 
  #try depth contours
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-100), size=c(0.3), colour="grey") +
  geom_contour(data = bath_cont_fort, aes(x=x, y=y, z=z), breaks=c(-50), size=c(0.3), colour="grey")+
  #ud50 for all whales combined in GOM
  geom_polygon(data = hwoct18,
               aes(x = long, y = lat, group = group),inherit.aes = FALSE,
               fill='red', alpha=0.25) +
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

months18telem <- ggarrange(jul18 + rremove("ylab") + rremove("xlab"),
                      aug18 + rremove("ylab") + rremove("xlab"),
                      sep18 + rremove("ylab") + rremove("xlab"),
                      oct18 + rremove("ylab") + rremove("xlab"),
                      
                      labels = c("A", "B", "C", "D"),
                      ncol=2,nrow=2) + 
  theme(plot.margin = margin(1,0.05,1,0.05, "cm"))
annotate_figure(months18telem, 
                left = text_grob("Latitude",rot=90, vjust=0.5),
                bottom = text_grob("Longitude", vjust=-4.0))
#can use annotate_figure 
ggsave("./figures/2018_monthly_overlap_telem.tiff",dpi=300,width=6,height=6,units="in")

