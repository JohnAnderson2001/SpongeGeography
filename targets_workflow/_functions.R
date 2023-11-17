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

infer_depths <- function(sponges_data_for_idw, reduced_data) {

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
		target_depth <- tidyph$transforms$depth$depth[(which.min(abs(tidyph$transforms$depth$depth - abs(sponges$idw_depths[sponge_index]))))]
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


extract_silica <- function(sponges, tidysilica) {

	sponges$silica <- NA

	for(sponge_index in sequence(nrow(sponges))) {
   	    #target_latitude <- tidyph$transforms$latitude$latitude[(which.min(abs(tidyph$transforms$latitude$latitude - sponges$decimalLatitude[sponge_index])))]
    	#target_longitude <- tidyph$transforms$longitude$longitude[(which.min(abs(tidyph$transforms$longitude$longitude - sponges$decimalLongitude[sponge_index])))]
    	target_depth <- tidysilica$transforms$depth$depth[(which.min(abs(tidysilica$transforms$depth$depth - abs(sponges$idw_depths[sponge_index]))))]
    	bound <- 0.25
    	sils <- tidysilica %>% hyper_filter(
        	latitude = between(latitude, sponges$decimalLatitude[sponge_index]-bound, sponges$decimalLatitude[sponge_index]+bound),
        	longitude = between(longitude, sponges$decimalLongitude[sponge_index]-bound, sponges$decimalLongitude[sponge_index]+bound),
        	depth = between(depth, target_depth-10, target_depth+10)
        	) %>% hyper_array()
   	    sponges$silica[sponge_index] <- median(sils$si, na.rm=TRUE)
    	cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(sils$si), " values")

	}
	return(sponges)
}


extract_temp <- function(sponges, tidytemp) {

	sponges$temperature <- NA

	for(sponge_index in sequence(nrow(sponges))) {
    	#target_latitude <- tidyph$transforms$latitude$latitude[(which.min(abs(tidyph$transforms$latitude$latitude - sponges$decimalLatitude[sponge_index])))]
    	#target_longitude <- tidyph$transforms$longitude$longitude[(which.min(abs(tidyph$transforms$longitude$longitude - sponges$decimalLongitude[sponge_index])))]
    	#target_depth <- tidysilica$transforms$depth$depth[(which.min(abs(tidysilica$transforms$depth$depth - sponges$depth[sponge_index])))]
    	bound <- 0.083
    	temps <- tidytemp %>% hyper_filter(
        	latitude = between(latitude, sponges$decimalLatitude[sponge_index]-bound, sponges$decimalLatitude[sponge_index]+bound),
        	longitude = between(longitude, sponges$decimalLongitude[sponge_index]-bound, sponges$decimalLongitude[sponge_index]+bound),
        	#depth = between(depth, target_depth-10, target_depth+10)
         	) %>% hyper_array()
    	sponges$temperature[sponge_index] <- median(temps$bottomT, na.rm=TRUE)
    	cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(temps$bottomT), " values")

	}
	return(sponges)
}

extract_viscosity <- function(sponges) {
	sponges$viscosity <- NA
	for(sponge_index in sequence(nrow(sponges))) {
		if((is.na(sponges$temperature[sponge_index]) == FALSE) && (sponges$temperature[sponge_index] >= 0.1)) {
			sponges$viscosity[sponge_index] <- ViscTD(Temp= sponges$temperature[sponge_index]+273.15, D= 1036)
		}
	}
	return(sponges)
}

#constructing models


datatree <- function(sponges, phy) {

	#prepare occurrence and environmental data; aggregate by genus
	cleandat <- sponges[,c("genus", "spicules", "ph", "temperature", "silica", "idw_depths", "viscosity"), ]
	phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
	tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
	sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
	depthdat <- aggregate(idw_depths ~ genus, FUN="mean", data=cleandat)
	depthdat$idw_depths <- abs(depthdat$idw_depths)
	viscdat <- aggregate(viscosity ~ genus, FUN="mean", data=cleandat)
	spicdat <- aggregate(spicules ~ genus, FUN="mean", data=cleandat)
	finalsponges <- merge(phdat, tempdat, by="genus")
	finalsponges <- merge(finalsponges, sildat, by="genus")
	finalsponges <- merge(finalsponges, depthdat, by="genus")
	finalsponges <- merge(finalsponges, viscdat, by="genus")
	finalsponges <- merge(finalsponges, spicdat, by="genus")
	row.names(finalsponges) <- finalsponges$genus

	#prepare phylogenetic tree to merge with occurrence/environmental data
	phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)
	RemoveDuplicateNames <- function(phy) {
		if(any(duplicated(phy$tip.label))) {
			phy$tip.label <- make.unique(phy$tip.label, sep = ".")
			cat("Duplicate names were found and have been renamed.\n")
		} else {
			cat("No duplicate names were found.\n")
		}
		return(phy)
	}

	prunedtree <- treedata(phy=RemoveDuplicateNames(phy), data=finalsponges)
	prunedtree$data <- as.data.frame(prunedtree$data)
	for (i in sequence(ncol(prunedtree$data))) {
		prunedtree$data[,i] <- as.numeric(prunedtree$data[,i])
	}
	return(prunedtree)
}


#use the original loops from Final_sponge_model_new_depth_data to extract all the environmental variables
#then use phyloglm for model because the data is binary
#use the same code for aggregating the mean conditions of each genus
#remember to  comment out the code reading the csv files because this causes problems when trying to run the code especially because the file path is so specific and will not be the same on every machine
#consider what it will mean if presence or absence does not appear influenced by any of the variables
#also make a map of where the spiculose and aspiculose sponges are located (should be similar to map from previous project)
#would boxplots be useful to show where most sponges fall on different environmental scales?

nobs.phyloglm <- function(object, ...) {
	return(object$n)
}

#possible function to compare models?
#compare_models <- function(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16) {
#	compare_1 <- rbind(ml$aic, m2$aic, m3$aic, m4$aic, m5$aic, m6$aic, m7$aic, m8$aic, m9$aic, m10$aic, m11$aic, m12$aic, m13$aic, m14$aic, m15$aic, m16$aic)
#	row.names(compare_1) <- c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16)
#	col.names(compare_!) <- "aic"
#	compare_1 %>% arrange(desc(aic))
#}

#binary data functions complete
#functions for building models with numeric complexity data:

create_complexity_data <- function(GBIF, complexity_dat) {
	reldat <- GBIF[,c("gbifID", "class", "order", "family", "genus", "species", "decimalLatitude", "decimalLongitude", "depth", "eventDate", "day", "month", "year")]
	colnames(complexity_dat)[1] = "species"
	prefinaldat <- merge(reldat, complexity_dat, by="species", all=T)
	finaldat <- prefinaldat[complete.cases(prefinaldat), ]
	finaldat$depth <- as.numeric(finaldat$depth)
	finaldat2 <- finaldat[!is.na(as.numeric(finaldat$decimalLatitude)), ]
	sponges <- finaldat2[!finaldat2$year < 1950, ]
	return(sponges)
}

complex_depth_extraction <- function(sponges, tidydepth) {
	sponges$newdepth <- NA
	sponges$decimalLatitude <- as.numeric(sponges$decimalLatitude)
	sponges$decimalLongitude <- as.numeric(sponges$decimalLongitude)

	sponges_data_for_idw <- data.frame()
	for(sponge_index in sequence(nrow(sponges))) {
		bound <- 0.00416667*10
			elevations <- tidydepth %>% hyper_filter(
				lat = between(lat, sponges$decimalLatitude[sponge_index]-bound, sponges$decimalLatitude[sponge_index]+bound),
				lon = between(lon, sponges$decimalLongitude[sponge_index]-bound, sponges$decimalLongitude[sponge_index]+bound),
				) %>% hyper_tibble()
				elevations$elevation[elevations$elevation>0] <- 0 # We need points on land so we can get shallow depths correct, but also don't want an island nearby to give a positive depth. So flatten the land.

				elevations_spatial <- as.data.frame(SpatialPointsDataFrame(elevations[,c('lon', 'lat')], elevations[,c('elevation')]))
				sponges_data_for_idw <- plyr::rbind.fill(sponges_data_for_idw, elevations_spatial)
				elevations_focal <- as.data.frame(SpatialPointsDataFrame(data.frame(lon=sponges$decimalLongitude[sponge_index], lat=sponges$decimalLatitude[sponge_index]), data=data.frame(elevation=NA)))
				sponges_data_for_idw <- plyr::rbind.fill(sponges_data_for_idw, elevations_focal)
				cat("\r", sponge_index, " of ", nrow(sponges), nrow(sponges_data_for_idw))
			
	}

	sponges_data_for_idw <- dplyr::distinct(sponges_data_for_idw)
	sponges_data_for_idw$gbif <- FALSE
	sponges_data_to_include <- sponges
	sponges_data_to_include$lon <- sponges_data_to_include$decimalLongitude
	sponges_data_to_include$lat <- sponges_data_to_include$decimalLatitude
	sponges_data_to_include$elevation <- sponges_data_to_include$depth
	sponges_data_to_include$gbif <- TRUE
	sponges_data_for_idw <- plyr::rbind.fill(sponges_data_for_idw, sponges_data_to_include)

	training_data <- subset(sponges_data_for_idw, gbif==FALSE & !is.na(elevation))
	test_data <- subset(sponges_data_for_idw, gbif==TRUE & !is.na(elevation))

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
	sponges$idw_depths <- new_depths
	sponges$idw_depths <- abs(sponges$idw_depths)
	return(sponges)
}

complex_ph_extraction <- function(sponges, tidyph) {
	sponges$ph <- NA

	for(sponge_index in sequence(nrow(sponges))) {
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

complex_silica_extraction <- function(sponges, tidysilica) {
	sponges$silica <- NA

	for(sponge_index in sequence(nrow(sponges))) {
		target_depth <- tidysilica$transforms$depth$depth[(which.min(abs(tidysilica$transforms$depth$depth - sponges$idw_depths[sponge_index])))]
		bound <- 0.25
		sils <- tidysilica %>% hyper_filter(
			latitude = between(latitude, sponges$decimalLatitude[sponge_index]-bound, sponges$decimalLatitude[sponge_index]+bound),
			longitude = between(longitude, sponges$decimalLongitude[sponge_index]-bound, sponges$decimalLongitude[sponge_index]+bound),
			depth = between(depth, target_depth-10, target_depth+10)
			) %>% hyper_array()
		sponges$silica[sponge_index] <- median(sils$si, na.rm=TRUE)
		cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(sils$si), " values")
	}
	return(sponges)
}

complex_temperature_extraction <- function(sponges, tidytemp) {
	sponges$temperature <- NA

	for(sponge_index in sequence(nrow(sponges))) {
		bound <- 0.083
		temps <- tidytemp %>% hyper_filter(
			latitude = between(latitude, sponges$decimalLatitude[sponge_index]-bound, sponges$decimalLatitude[sponge_index]+bound),
			longitude = between(longitude, sponges$decimalLongitude[sponge_index]-bound, sponges$decimalLongitude[sponge_index]+bound),
			) %>% hyper_array()
		sponges$temperature[sponge_index] <- median(temps$bottomT, na.rm=TRUE)
		cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(temps$bottomT), " values")
	}
	return(sponges)
}

complex_viscosity_extraction <- function(sponges) {
	sponges$viscosity <- NA
	for(sponge_index in sequence(nrow(sponges))) {
		if((is.na(sponges$temperature[sponge_index]) == FALSE) && (sponges$temperature[sponge_index] >= 0.1)) {
			sponges$viscosity[sponge_index] <- ViscTD(Temp= sponges$temperature[sponge_index]+273.15, D= 1036)
		}
	}
	return(sponges)
}

complex_datatree <- function(sponges, phy) {
	cleandat <- sponges[,c("genus", "SpiculeTypes", "ph", "temperature", "silica", "idw_depths", "viscosity"), ]
	phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
	tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
	sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
	depthdat <- aggregate(idw_depths ~ genus, FUN="mean", data=cleandat)
	depthdat$idw_depths <- abs(depthdat$idw_depths)
	viscdat <- aggregate(viscosity ~ genus, FUN="mean", data=cleandat)
	spicdat <- aggregate(SpiculeTypes ~ genus, FUN="mean", data=cleandat)
	finalsponges <- merge(phdat, tempdat, by="genus")
	finalsponges <- merge(finalsponges, sildat, by="genus")
	finalsponges <- merge(finalsponges, depthdat, by="genus")
	finalsponges <- merge(finalsponges, viscdat, by="genus")
	finalsponges <- merge(finalsponges, spicdat, by="genus")
	row.names(finalsponges) <- finalsponges$genus

	phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)

	RemoveDuplicateNames <- function(phy) {
		if(any(duplicated(phy$tip.label))) {
			phy$tip.label <- make.unique(phy$tip.label, sep = ".")
			cat("Duplicate names were found and have been renamed.\n")
		} else {
			cat("No duplicate names were found.\n")
		}
		return(phy)
	}

	prunedtree <- treedata(phy=RemoveDuplicateNames(phy), data=finalsponges)

	prunedtree$data <- as.data.frame(prunedtree$data)
	for (i in sequence(ncol(prunedtree$data))) {
		prunedtree$data[,i] <- as.numeric(prunedtree$data[,i])
	}
	prunedtree$data$LogSpiculeTypes <- log(prunedtree$data$SpiculeTypes)

	return(prunedtree)
}

complex_bm_model <- function(prunedtree) {
	phylomodelbmgenus <- phylolm(LogSpiculeTypes ~ ph + temperature + silica + idw_depths + viscosity, data=prunedtree$data, phy=prunedtree$phy, model="BM")
	bestmodelbmgenus <- dredge(phylomodelbmgenus)
	return(bestmodelbmgenus)
}

complex_randomroot_model <- function(prunedtree) {
	phylomodelrandomrootgenus <- phylolm(LogSpiculeTypes ~ ph + temperature + silica + idw_depths + viscosity, data=prunedtree$data, phy=prunedtree$phy, model="OUrandomRoot")
	bestmodelrandomrootgenus <- dredge(phylomodelrandomrootgenus)
	return(bestmodelrandomrootgenus)
}

#make adjusted models with different predictors and transformations on data
complex_data_adjustments <- function(spongedata) {
	#divide sponges with silica and calcite spicules so sponges that do not use silica can be accounted for in the model
	#calcarea <- sponges[sponges$class %in% "Calcarea", ]
	#notcalcarea <- sponges[!sponges$class %in% "Calcarea", ]
	#calcarea$silicaspicules <- 0
	#notcalcarea$silicaspicules <- 1
	#adjusted_data <- rbind(notcalcarea, calcarea)
	sponges <- spongedata
	sponges$silicaspicules <- 1
	for (row in sequence(nrow(sponges))) {
		if (sponges$class[row] == "Calcarea") {
			sponges$silicaspicules[row] <- 0
		}
	}
	#divide photic and non-photic depth zones
	#find source for photic zone being 200m
	sponges$photic <- 0
	for (row in sequence(nrow(sponges))) {
		if (sponges$idw_depths[row] < 200) {
			sponges$photic[row] <- 1
		}
	}
	sponges$nonphotic <- 0
	for (row in sequence(nrow(sponges))) {
		if (sponges$idw_depths[row] > 200) {
			sponges$nonphotic[row] <- 1
		}
	}

	return(sponges)
}

complex_tree_adjustments <- function(sponges, phy) {
	cleandat <- sponges[,c("genus", "SpiculeTypes", "ph", "temperature", "silica", "idw_depths", "viscosity", "silicaspicules", "photic", "nonphotic"), ]
	phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
	tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
	sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
	depthdat <- aggregate(idw_depths ~ genus, FUN="mean", data=cleandat)
	viscdat <- aggregate(viscosity ~ genus, FUN="max", data=cleandat)
	silicaspicsdat <- aggregate(silicaspicules ~ genus, FUN="max", data=cleandat)
	photicdat <- aggregate(photic ~ genus, FUN="max", data=cleandat)
	nonphoticdat <- aggregate(nonphotic ~ genus, FUN="max", data=cleandat)
	depthdat$idw_depths <- abs(depthdat$idw_depths)
	spicdat <- aggregate(SpiculeTypes ~ genus, FUN="mean", data=cleandat)
	finalsponges <- merge(phdat, tempdat, by="genus")
	finalsponges <- merge(finalsponges, sildat, by="genus")
	finalsponges <- merge(finalsponges, depthdat, by="genus")
	finalsponges <- merge(finalsponges, viscdat, by="genus")
	finalsponges <- merge(finalsponges, silicaspicsdat, by="genus")
	finalsponges <- merge(finalsponges, photicdat, by="genus")
	finalsponges <- merge(finalsponges, nonphoticdat, by="genus")
	finalsponges <- merge(finalsponges, spicdat, by="genus")
	row.names(finalsponges) <- finalsponges$genus

	phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)

	RemoveDuplicateNames <- function(phy) {
		if(any(duplicated(phy$tip.label))) {
			phy$tip.label <- make.unique(phy$tip.label, sep = ".")
			cat("Duplicate names were found and have been renamed.\n")
		} else {
			cat("No duplicate names were found.\n")
		}
		return(phy)
	}

	prunedtree <- treedata(phy=RemoveDuplicateNames(phy), data=finalsponges)

	prunedtree$data <- as.data.frame(prunedtree$data)
	for (i in sequence(ncol(prunedtree$data))) {
		prunedtree$data[,i] <- as.numeric(prunedtree$data[,i])
	}
	prunedtree$data$LogSpiculeTypes <- log(prunedtree$data$SpiculeTypes)

	return(prunedtree)
}

complex_bm_model_adjustment <- function(prunedtree) {
	phylomodelbmgenus <- phylolm(LogSpiculeTypes ~ ph + temperature + silica + log(idw_depths + 1) + viscosity + photic + nonphotic * silicaspicules, data=prunedtree$data, phy=prunedtree$phy, model="BM")
	bestmodelbmgenus <- dredge(phylomodelbmgenus)
	return(bestmodelbmgenus)
	#is the log depth transformation correct?
}

complex_rr_model_adjustment <- function(prunedtree) {
	phylomodelrrgenus <- phylolm(LogSpiculeTypes ~ ph + temperature + silica + log(idw_depths + 1) + viscosity + photic + nonphotic * silicaspicules, data=prunedtree$data, phy=prunedtree$phy, model="OUrandomRoot")
	bestmodelrrgenus <- dredge(phylomodelrrgenus)
	return(bestmodelrrgenus)
	#is the log depth transformation correct?
}

binary_data_adjustments <- function(spongedata) {
	#divide sponges with silica and calcite spicules so sponges that do not use silica can be accounted for in the model
	#calcarea <- sponges[sponges$class %in% "Calcarea", ]
	#notcalcarea <- sponges[!sponges$class %in% "Calcarea", ]
	#calcarea$silicaspicules <- 0
	#notcalcarea$silicaspicules <- 1
	#adjusted_data <- rbind(notcalcarea, calcarea)
	sponges <- spongedata
	sponges$idw_depths <- abs(sponges$idw_depths)
	sponges$silicaspicules <- 1
	for (row in sequence(nrow(sponges))) {
		if (sponges$class[row] == "Calcarea") {
			sponges$silicaspicules[row] <- 0
		}
		if (sponges$spicules[row] == 0) {
			sponges$silicaspicules[row] <- 0
		}
	}
	#divide photic and non-photic depth zones
	#find source for photic zone being 200m
	sponges$photic <- 0
	for (row in sequence(nrow(sponges))) {
		if (sponges$idw_depths[row] < 200) {
			sponges$photic[row] <- 1
		}
	}
	sponges$nonphotic <- 0
	for (row in sequence(nrow(sponges))) {
		if (sponges$idw_depths[row] > 200) {
			sponges$nonphotic[row] <- 1
		}
	}

	return(sponges)
}

#functions for sister-group analysis for presence/absence of spicules
sisters_analysis <- function(sponges, phy, cutoff=2.1) {
	#cleaned <- sis_clean(phy=RemoveDuplicateNames(phy), complexity_depths_ph_silica_and_temp, first_col_names = TRUE)
	#use pruned tree instead?
	cleandat <- sponges[,c("genus", "SpiculeTypes", "ph", "temperature", "silica", "idw_depths", "viscosity"), ]
	phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
	tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
	sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
	depthdat <- aggregate(idw_depths ~ genus, FUN="mean", data=cleandat)
	depthdat$idw_depths <- abs(depthdat$idw_depths)
	viscdat <- aggregate(viscosity ~ genus, FUN="mean", data=cleandat)
	spicdat <- aggregate(SpiculeTypes ~ genus, FUN="mean", data=cleandat)
	finalsponges <- merge(phdat, tempdat, by="genus")
	finalsponges <- merge(finalsponges, sildat, by="genus")
	finalsponges <- merge(finalsponges, depthdat, by="genus")
	finalsponges <- merge(finalsponges, viscdat, by="genus")
	finalsponges <- merge(finalsponges, spicdat, by="genus")
	#inalsponges$species <- sub(" ", "_", finalsponges$species)
	row.names(finalsponges) <- finalsponges$genus

	phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)
	RemoveDuplicateNames <- function(phy) {
		if(any(duplicated(phy$tip.label))) {
			phy$tip.label <- make.unique(phy$tip.label, sep = ".")
			cat("Duplicate names were found and have been renamed.\n")
		} else {
			cat("No duplicate names were found.\n")
		}
		return(phy)
	}

	cleaned <- sis_clean(RemoveDuplicateNames(phy), finalsponges)
	#print(colnames(cleaned$traits))
	phy <- cleaned$phy
	traits <- cleaned$traits
	print(quantile(as.numeric(traits[,"SpiculeTypes"])))
	trait <- sis_discretize(as.numeric(traits[,"SpiculeTypes"]), cutoff=cutoff, use_percentile=FALSE)
	names(trait) <- rownames(traits)
	#print(quantile(trait))
	#trait <- cleaned$traits[,6]
	sisters <- sis_get_sisters(phy)
	sisters_comparison <- sis_format_comparison(sisters, trait, phy)
	concatenate_states <- function(x) {
		return(paste0(sort(unique(unname(x[[1]]))), collapse=""))
	}
	sisters_comparison$left.trait.simplified <- NA
	sisters_comparison$right.trait.simplified <- NA
	for (i in sequence(nrow(sisters_comparison))) {
		sisters_comparison$left.trait.simplified[i] <- concatenate_states(sisters_comparison$left.trait[i])
		sisters_comparison$right.trait.simplified[i] <- concatenate_states(sisters_comparison$right.trait[i])
	}

	sisters_comparison_informative <- subset(sisters_comparison, left.trait.simplified!=right.trait.simplified & nchar(left.trait.simplified)==1 & nchar(right.trait.simplified)==1)
	two_taxon_pairs <- data.frame()
	all_pairs <- data.frame()
	for (sister_index in sequence(nrow(sisters_comparison_informative))) {
		if(sisters_comparison_informative$ntax.total[sister_index]==2) {
			focal_left <- phy$tip.label[unlist(sisters_comparison_informative$left[sister_index])]
			focal_right <- phy$tip.label[unlist(sisters_comparison_informative$right[sister_index])]
			paired_traits <-as.data.frame(traits[rownames(traits) %in% c(focal_left, focal_right),])
			paired_traits$pair_node <- sisters_comparison_informative$node[sister_index]
			two_taxon_pairs <- dplyr::bind_rows(two_taxon_pairs, paired_traits)
			all_pairs <- dplyr::bind_rows(all_pairs, paired_traits)
		} else {
			focal_left <- phy$tip.label[unlist(sisters_comparison_informative$left[sister_index])]	
			focal_right <- phy$tip.label[unlist(sisters_comparison_informative$right[sister_index])]
			focal_left_traits <- as.data.frame(traits[rownames(traits) %in% c(focal_left),])
			focal_right_traits <- as.data.frame(traits[rownames(traits) %in% c(focal_right),])
			if(ncol(focal_left_traits)==1) {
				focal_left_traits <- as.data.frame(t(focal_left_traits))
				rownames(focal_left_traits) <- NULL
			}
			if(ncol(focal_right_traits)==1) {
				focal_right_traits <- as.data.frame(t(focal_right_traits))
				rownames(focal_right_traits) <- NULL
			}
			focal_left_concat <- data.frame(matrix(nrow=1,ncol=ncol(focal_left_traits)))
			focal_right_concat <- data.frame(matrix(nrow=1,ncol=ncol(focal_right_traits)))
			for(col_index in sequence(ncol(focal_left_concat))) {
				focal_left_concat[,col_index] <- paste0(sort(unique(unname(focal_left_traits[,col_index]))), collapse=", ")	
				focal_right_concat[,col_index] <- paste0(sort(unique(unname(focal_right_traits[,col_index]))), collapse=", ")
			}
			paired_traits <- dplyr::bind_rows(focal_left_concat, focal_right_concat)
			colnames(paired_traits) <- colnames(focal_left_traits)
			paired_traits$pair_node <- sisters_comparison_informative$node[sister_index]
			all_pairs <- dplyr::bind_rows(all_pairs, paired_traits)
		}
	}
	original_width <- ncol(all_pairs)
	for (i in 2:original_width){
		all_pairs[,ncol(all_pairs)+1] <- rep(NA, nrow(all_pairs))
		all_pairs[,ncol(all_pairs)+1] <- rep(NA, nrow(all_pairs))
		for (row_index in sequence(nrow(all_pairs))) {
			values <- as.numeric(strsplit(as.character(all_pairs[row_index,i]), ", ")[[1]])
			all_pairs[row_index,-1+ncol(all_pairs)] <- min(values)
			all_pairs[row_index,ncol(all_pairs)] <- max(values)
		}
		colnames(all_pairs)[ncol(all_pairs)-1] <- paste0(colnames(all_pairs)[i], "_min")
		colnames(all_pairs)[ncol(all_pairs)] <- paste0(colnames(all_pairs)[i], "_max")
	}
	rownames(all_pairs) <- NULL
	all_pairs <- all_pairs[order(all_pairs$pair_node, all_pairs$SpiculeTypes_min),]
	return(all_pairs)
	#pairs <- sis_format_simpified(sisters_comparison)
	#return(pairs)
	#sis_test(pairs)
	#test <- sis_get_sisters(phy)
	#test2 <- sis_format_comparison(sisters=test, trait=pruned_tree$data, phy=phy)
}

do_many_cutoffs <- function(sponges, phy) {
	results <- data.frame()
	cutoffs <- seq(from=0.1, to=14, by=0.1)
	for (cutoff in cutoffs) {
		try({
		sister_results <- sisters_analysis(sponges, phy, cutoff=cutoff)
		sister_results$cutoff <- cutoff
		results <- dplyr::bind_rows(results, sister_results)
		})
	}	
	return(results)
}

summarize_cutoffs <- function(complex_sister_analysis) {
	results <- complex_sister_analysis |> group_by(cutoff) |> summarize(npairs =n()/2)	
	return(results)
}

binomial_complexity_prep <- function(sponges, phy) {
	for(i in sequence(nrow(sponges))) {
		if (sponges$SpiculeTypes[i] >= 2.1) {
			sponges$SpiculeTypes[i] <- 1
		} else {
			sponges$SpiculeTypes[i] <- 0
		}
	}

	cleandat <- sponges[,c("genus", "SpiculeTypes", "ph", "temperature", "silica", "idw_depths", "viscosity"), ]
	phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
	tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
	sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
	depthdat <- aggregate(idw_depths ~ genus, FUN="mean", data=cleandat)
	depthdat$idw_depths <- abs(depthdat$idw_depths)
	viscdat <- aggregate(viscosity ~ genus, FUN="mean", data=cleandat)
	spicdat <- aggregate(SpiculeTypes ~ genus, FUN="mean", data=cleandat)
	finalsponges <- merge(phdat, tempdat, by="genus")
	finalsponges <- merge(finalsponges, sildat, by="genus")
	finalsponges <- merge(finalsponges, depthdat, by="genus")
	finalsponges <- merge(finalsponges, viscdat, by="genus")
	finalsponges <- merge(finalsponges, spicdat, by="genus")
	#inalsponges$species <- sub(" ", "_", finalsponges$species)
	row.names(finalsponges) <- finalsponges$genus

	phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)
	RemoveDuplicateNames <- function(phy) {
		if(any(duplicated(phy$tip.label))) {
			phy$tip.label <- make.unique(phy$tip.label, sep = ".")
			cat("Duplicate names were found and have been renamed.\n")
		} else {
			cat("No duplicate names were found.\n")
		}
		return(phy)
	}

	prunedtree <- treedata(phy=RemoveDuplicateNames(phy), data=finalsponges)

	prunedtree$data <- as.data.frame(prunedtree$data)
	for (i in sequence(ncol(prunedtree$data))) {
		prunedtree$data[,i] <- as.numeric(prunedtree$data[,i])
	}
	#prunedtree$data$LogSpiculeTypes <- log(prunedtree$data$SpiculeTypes)
	prunedtree$data = prunedtree$data[(prunedtree$data$SpiculeTypes == 1 | prunedtree$data$SpiculeTypes == 0), ]
	#for(i in sequence(nrow(prunedtree$data))) {
		#if(prunedtree$data$SpiculeTypes[i] != 1 | prunedtree$data$SpiculeTypes[i] != 0) {
			#prunedtree$data <- prunedtree$data[!prunedtree$data[i],]
		#}
	#}

	return(prunedtree)
}


do_all_pics <- function(sponges, phy) {
	cleandat <- sponges[,c("genus", "SpiculeTypes", "ph", "temperature", "silica", "idw_depths", "viscosity"), ]
	phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
	tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
	sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
	depthdat <- aggregate(idw_depths ~ genus, FUN="mean", data=cleandat)
	depthdat$idw_depths <- abs(depthdat$idw_depths)
	spicdat <- aggregate(SpiculeTypes ~ genus, FUN="mean", data=cleandat)
	viscositydat <- aggregate(viscosity ~ genus, FUN="mean", data=cleandat)
	finalsponges <- merge(phdat, tempdat, by="genus")
	finalsponges <- merge(finalsponges, sildat, by="genus")
	finalsponges <- merge(finalsponges, depthdat, by="genus")
	finalsponges <- merge(finalsponges, spicdat, by="genus")
	finalsponges <- merge(finalsponges, viscositydat, by="genus")
	#inalsponges$species <- sub(" ", "_", finalsponges$species)
	row.names(finalsponges) <- finalsponges$genus

	phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)
	RemoveDuplicateNames <- function(phy) {
		if(any(duplicated(phy$tip.label))) {
			phy$tip.label <- make.unique(phy$tip.label, sep = ".")
			cat("Duplicate names were found and have been renamed.\n")
		} else {
			cat("No duplicate names were found.\n")
		}
		return(phy)
	}

	cleaned <- sis_clean(RemoveDuplicateNames(phy), finalsponges)
	#print(colnames(cleaned$traits))
	phy <- cleaned$phy
	traits <- cleaned$traits
	traits <- traits[,!grepl('genus', colnames(traits))]
	traits <- as.data.frame(traits)
	traits$logSpiculeTypes <- log(as.numeric(traits$SpiculeTypes))
	results <- data.frame()
	for (i in sequence(ncol(traits))) {
		comparison <- t(t(ape::pic(as.numeric(traits[,i]), phy)))
		if(i==1) {
			results <- comparison
		} else {
			results <- cbind(results, comparison)
		}
		colnames(results)[ncol(results)] <- colnames(traits)[i]	
	}
	return(as.data.frame(results))
}

compute_pics_summaries <- function(pics, focal="logSpiculeTypes") {
	options(na.action = "na.fail")
	formula <- paste0(focal, " ~ ph + temperature + silica + idw_depths + viscosity + 0")
	# lm regression
	lm_results <- lm(formula, data=pics)
	return(MuMIn::dredge(lm_results))
}

#find sponges in GBIF
#relevanttest <- complexityGBIF[,c("gbifID", "class", "order", "family", "genus", "species", "decimalLatitude", "decimalLongitude", "depth", "year")]
#test <- relevanttest[(complexityGBIF$species %in% "Galaxia gaviotensis"), ] #replace string with the name of a species to see if it is in GBIF

#testnotree <- glm(spicules ~ idw_depths + ph + silica + temperature, data=depths_pH_silica_and_temp, family="binomial")
#step(testnotree)

