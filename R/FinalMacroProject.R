test <- read.delim("0196959-210914110416597.csv")

cleandat <- test [ !test$year %in% NA, ]
cleandat2 <- cleandat [ !cleandat$individualCount %in% NA, ]
cleandat3 <- cleandat2 [ cleandat2$year %in% 1993, ]

cleandat2$ones <- rep(1, nrow(cleandat2))

ag <- aggregate(ones ~ year, FUN="sum", data=cleandat2)


test2 <- read.delim("0197294-210914110416597.csv")

cleandat <- test2 [ !test2$year %in% NA, ]
cleandat2 <- test2 [ !cleandat$individualCount %in% NA, ]

#use occurrence data from a range of time so long as it overlaps / is close to the range of time for the environmental data
#for the first analysis, just find the correlation between environmental variables and the occurrence of spiculeless sponges and sponges with spicules
#find average environmental conditions for species with spicules and average conditions for species without to find out what kind of environments they thrive in
#can do a temporal analysis as the second step if the data is available


#ocurdat <- read.delim("0224427-210914110416597.csv")
#takes long time to load, run this in regular R and not RStudio
#don't use R integration for Bio-ORACLE because it doesn't work

demos <- ocurdat [ocurdat$class %in% 'Demospongiae', ]
demoscleaned <- demos[!is.na(demos$species),]

cleandepthdat <- ocurdat[!ocurdat$depth %in% NA, ]
reldat <- cleandepthdat[,c("gbifID", "class", "order", "family", "genus", "species", "decimalLatitude", "decimalLongitude", "depth")]
demos <- reldat[reldat$class %in% 'Demospongiae', ]
nospicules <- demos[demos$order %in% c('Dendroceratida', 'Dictyoceratida', 'Verongiida', 'Chondrosiida'), ]
finalnospics <- nospicules[!nospicules$genus %in% 'Chondrilla', ]

finalspics <- read.csv("finalspics.csv")
finalnospics <-read.csv("finalnospics.csv")

pH <- raster("Present.Surface.pH.BOv2_2.asc")
silica <- raster("Present.Surface.Silicate.Mean.asc")

#the extremely important extraction:
points(finalnospics$decimalLongitude, finalnospics$decimalLatitude, col="darkblue", pch=19)

finalnospics$pH <- extract(pH, finalnospics[,c("decimalLongitude","decimalLatitude")], method="simple")



depths <- c(finalspics$depth, finalnospics$depth)
max(depths)
mean(depths)
q25 <- quantile(depths)[2]

test <- depths[depths < q25]

#TAXA WITH NO SPICULES:
#Dendroceratida, Dictyoceratida, Verongiida, Chondrosiida (except Chondrilla)

#spicules are possessed by subclasses Haploscleromorpha and Heteroscleromorpha


#a test of trying to get data into raster
nospics$presence <- rep(1, nrow(nospics))
testdat <- nospics[,c("decimalLatitude", "decimalLongitude", "presence")] #doesn't work



# establish rough extent of study sites						
ROI <- extent(-130,-50,-10,50)
# plot temperature grids over this extent for the year 1959
plot(crop(tempbrick1.rotated[[1]],ROI))						
# add points of sites
points(latlon.df$Lon, latlon.df$Lat, col="darkblue", pch=19)	## looks right!

# now extract the points from first year of dataset 1959
temp.df <- extract(tempbrick1.rotated, latlon.df[,c("Lon","Lat")], method="simple")	# note NAs for palo and guanica

