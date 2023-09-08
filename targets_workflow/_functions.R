create_reduced_data <- function(tidydepth) {
	if(file.exists(file.path("data", "reduced_data.csv"))) {
		reduced_data <- read.csv(file.path("data", "reduced_data.csv"))
		reduced_data$decimalLatitude <- as.numeric(reduced_data$decimalLatitude)
		reduced_data$decimalLongitude <- as.numeric(reduced_data$decimalLongitude)

		return(reduced_data)
	} else {
		ocurdat <- read.delim("C:\\Users\\jjera\\Documents\\RStuff\\0163120-230224095556074.csv")
		head(ocurdat)
		reldat <- ocurdat[,c("gbifID", "class", "order", "family", "genus", "species", "decimalLatitude", "decimalLongitude", "depth", "eventDate", "day", "month", "year")]
		nospicules <- reldat[reldat$order %in% c('Dendroceratida', 'Dictyoceratida', 'Verongiida', 'Chondrosiida'), ]
		finalnospics <- nospicules[!nospicules$genus %in% 'Chondrilla', ]
		finalspics <- reldat[!reldat$order %in% c('Dendroceratida', 'Dictyoceratida', 'Verongiida', 'Chondrosiida'), ]

		finalnospics$spicules <- 0
		finalspics$spicules <- 1

		finaldat <- rbind(finalnospics, finalspics) #731777 rows
		#write.csv(finaldat, "finalbinaryspongedata.csv")

		dat <- finaldat [complete.cases(finaldat), ] #557357 rows
		dat2 <- dat[!is.na(as.numeric(dat$decimalLatitude)), ] #540327 rows; gets rid of bad location data
		dat3 <- dat2[!is.na(as.numeric(dat2$decimalLongitude)), ] #540327 rows
		sponges <- dat3[!dat3$year < 1950, ] #531773 rows
		sponges <- sponges[nchar(sponges$genus)>2,]
		#write.csv(sponges, "finalbinaryspongedata.csv")

		#sponges <- read.csv("C:\\Users\\jjera\\Documents\\RStuff\\finalbinaryspongedata.csv")

		#extracting depth


		sponges$newdepth <- NA
		sponges$npoints_newdepth <- NA

		sponges_data_for_idw <- data.frame()
		reduced_data <- data.frame()

		all_sponge_genera <- unique(sponges$genus)
		#new code to get all points

		#take a subset of the data; no more than 10 of each species within each genus (and no more than 10 unnamed species within each genus)
		for(sponge_genus_index in sequence(length(all_sponge_genera))) {
			sponge_genus <- all_sponge_genera[sponge_genus_index]
			genus_subset <- subset(sponges, genus==sponge_genus)
			for (sponge_species in unique(genus_subset$species)) {
				species_subset <- subset(genus_subset, species==sponge_species)
				if(nrow(species_subset)>10) {
					species_subset <- dplyr::sample_n(species_subset, 10)
				}
					reduced_data <- plyr::rbind.fill(reduced_data, species_subset)
			}
			cat("\r", sponge_genus_index, " of ", length(all_sponge_genera))
		}
		reduced_data$decimalLatitude <- as.numeric(reduced_data$decimalLatitude)
		reduced_data$decimalLongitude <- as.numeric(reduced_data$decimalLongitude)

		return(reduced_data)
	}
}

extracts_depths_around_sponges <- function(reduced_data, tidydepth, bound=0.00416667*5) {
	sponges_data_for_idw <- data.frame()
#depth
	for(sponge_index in sequence(nrow(reduced_data))) {
		elevations <- tidydepth %>% hyper_filter(
			lat = between(lat, as.numeric(reduced_data$decimalLatitude[sponge_index])-bound, as.numeric(reduced_data$decimalLatitude[sponge_index])+bound),
			lon = between(lon, as.numeric(reduced_data$decimalLongitude[sponge_index])-bound, as.numeric(reduced_data$decimalLongitude[sponge_index])+bound),
		#depth = between(depth, target_depth-10, target_depth+10)
		) %>% hyper_tibble()
		elevations$elevation[elevations$elevation>0] <- 0 # We need points on land so we can get shallow depths correct, but also don't want an island nearby to give a positive depth. So flatten the land.

		elevations_spatial <- as.data.frame(SpatialPointsDataFrame(elevations[,c('lon', 'lat')], elevations[,c('elevation')]))
		sponges_data_for_idw <- plyr::rbind.fill(sponges_data_for_idw, elevations_spatial)
		elevations_focal <- as.data.frame(SpatialPointsDataFrame(data.frame(lon=reduced_data$decimalLongitude[sponge_index], lat=reduced_data$decimalLatitude[sponge_index]), data=data.frame(elevation=NA)))
		sponges_data_for_idw <- plyr::rbind.fill(sponges_data_for_idw, elevations_focal)
		cat("\r", sponge_index, " of ", nrow(reduced_data))
	}
	return(sponges_data_for_idw)
}

infer_depths <- function(sponges_data_for_idw) {

	sponges_data_for_idw <- dplyr::distinct(sponges_data_for_idw)
	sponges_data_for_idw$gbif <- FALSE
	#sponges_data_for_idw$lon <- as.numeric(sponges_data_for_idw$lon)
	#sponges_data_for_idw$lat <- as.numeric(sponges_data_for_idw$lat)
	sponges_data_to_include <- reduced_data
	sponges_data_to_include$lon <- sponges_data_to_include$decimalLongitude
	sponges_data_to_include$lat <- sponges_data_to_include$decimalLatitude
	sponges_data_to_include$elevation <- sponges_data_to_include$depth
	#sponges_data_to_include$depth <- as.numeric(sponges_data_to_include$depth)
	sponges_data_to_include$gbif <- TRUE
	sponges_data_for_idw <- plyr::rbind.fill(sponges_data_for_idw, sponges_data_to_include)

	training_data <- subset(sponges_data_for_idw, gbif==FALSE & !is.na(elevation))
	training_data$elevation <- as.numeric(training_data$elevation)
	#training_data$elevation <- !is.na(training_data$elevation)
	test_data <- subset(sponges_data_for_idw, gbif==TRUE & !is.na(elevation))
	test_data$elevation <- as.numeric(test_data$elevation)
	#test_data$elevation <- !is.na(test_data$elevation)

	f1 <- function(x, test, train) {
		nmx <- x[1]
		idp <- x[2]
		if(nmx<1 | idp < 0.001) {return(Inf)}
		m <- gstat(formula=elevation~1, locations=~lon+lat, data=train, nmax=nmx, set=list(idp=idp))
		p <- predict(m, newdata=test, depug.level=0)$var1.pred
		return(Metrics::rmse(test$elevation, p))
	} 

	opt <- optim(c(8, 0.5), f1, test=test_data, train=training_data)

	m <- gstat(formula=elevation~1, locations=~lon+lat, data=subset(sponges_data_for_idw, !is.na(elevation) & gbif==FALSE), nmax=opt$par[1], set=list(idp=opt$par[2]))
	new_depths <- predict(m, newdata=subset(sponges_data_for_idw, gbif==TRUE))$var1.pred
	#start here
	reduced_data$idw_depths <- new_depths
	#plot(log1p(abs(sponges$depth)), log1p(abs(sponges$idw_depths)))
	#write.csv(sponges, "all_sponges_with_depths.csv")
	#code for optimization comes from https://rspatial.org/analysis/4-interpolation.html
	
	return(reduced_data)
	
}

extract_ph <- function(sponges, tidyph) {

	#extracting pH
	#tidyph <- tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\pH\\allpH.nc")
	sponges$ph <- NA

	for(sponge_index in sequence(nrow(sponges))) {
		#target_latitude <- tidyph$transforms$latitude$latitude[(which.min(abs(tidyph$transforms$latitude$latitude - sponges$decimalLatitude[sponge_index])))]
		#target_longitude <- tidyph$transforms$longitude$longitude[(which.min(abs(tidyph$transforms$longitude$longitude - sponges$decimalLongitude[sponge_index])))]
		target_depth <- tidyph$transforms$depth$depth[(which.min(abs(tidyph$transforms$depth$depth - sponges$idw_depths[sponge_index])))]
		bound <- 0.25
		phs <- tidyph %>% hyper_filter(
			latitude = between(latitude, sponges$decimalLatitude[sponge_index]-bound, sponges$decimalLatitude[sponge_index]+bound),
			longitude = between(longitude, sponges$decimalLongitude[sponge_index]-bound, sponges$decimalLongitude[sponge_index]+bound),
			depth = between(depth, target_depth-10, target_depth+10)
			) %>% hyper_array()
		sponges$ph[sponge_index] <- median(phs$ph, na.rm=TRUE)
		cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(phs$ph), " values")

	}
	return(sponges)
}





#use the original loops from Final_sponge_model_new_depth_data to extract all the environmental variables
#then use phyloglm for model because the data is binary
#use the same code for aggregating the mean conditions of each genus
#remember to  comment out the code reading the csv files because this causes problems when trying to run the code especially because the file path is so specific and will not be the same on every machine
#consider what it will mean if presence or absence does not appear influenced by any of the variables
#also make a map of where the spiculose and aspiculose sponges are located (should be similar to map from previous project)
#would boxplots be useful to show where most sponges fall on different environmental scales?
